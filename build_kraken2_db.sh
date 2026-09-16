#!/usr/bin/env bash
# Build a Kraken2 database from RefSeq for ancient-DNA taxonomic classification.
#
# Steps: taxonomy -> UniVec_Core -> download RefSeq genomes -> dustmask + tag ->
# add to library -> build the index. Downloading and masking run under GNU parallel.
#
# All settings come from the command line; see --help.

set -euo pipefail

SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- defaults ---------------------------------------------------------------
DB=$PWD/kraken2_refseq_k25
FASTAS=$PWD/fastas
WORK=$PWD/work
KMER=25
MINIMIZER=22
SPACES=0
MAXDB=1100000000000                     # 1.1 TB - must fit in RAM at classify time
THREADS=${SLURM_CPUS_PER_TASK:-8}
DL_JOBS=8                               # NCBI asks for <= 10 parallel connections
LIBS="archaea bacteria fungi invertebrate plant protozoa vertebrate_mammalian"
VIRAL=1                                 # all viral assemblies (includes PhiX174)
PLASMID=1                               # refseq/release/plasmid - needs the taxid maps
SKIP_MAPS=0                             # 1 = no accession2taxid download (breaks plasmids)
TEST_N=0
STORAGE=""

usage() {
cat <<EOF
usage: $(basename "$0") [options]

  --db DIR              database directory            (default: $DB)
  --fastas DIR          where genomes are downloaded  (default: $FASTAS)
  --kmer N              k-mer length                  (default: $KMER)
  --minimizer N         minimizer length, keep = kmer (default: $MINIMIZER)
  --spaces N            masked positions, keep 0      (default: $SPACES)
  --max-db-size BYTES   hash table cap                (default: $MAXDB)
  --threads N           masking + build threads       (default: $THREADS)
  --groups "a b c"      RefSeq groups to download     (default: $LIBS)
  --no-viral            skip the viral assemblies
  --no-plasmid          skip refseq/release/plasmid
  --skip-maps           do not download accession2taxid (only safe without plasmids)
  --test N              only the first N genomes, to check the recipe
  --storage DIR         rsync the finished database here
EOF
}

while [ $# -gt 0 ]; do
  case "$1" in
    --db)           DB=$2; shift 2;;
    --fastas)       FASTAS=$2; shift 2;;
    --kmer)         KMER=$2; shift 2;;
    --minimizer)    MINIMIZER=$2; shift 2;;
    --spaces)       SPACES=$2; shift 2;;
    --max-db-size)  MAXDB=$2; shift 2;;
    --threads)      THREADS=$2; shift 2;;
    --groups)       LIBS=$2; shift 2;;
    --no-viral)     VIRAL=0; shift;;
    --no-plasmid)   PLASMID=0; shift;;
    --skip-maps)    SKIP_MAPS=1; shift;;
    --test)         TEST_N=$2; shift 2;;
    --storage)      STORAGE=$2; shift 2;;
    -h|--help)      usage; exit 0;;
    *) echo "unknown option: $1" >&2; usage; exit 1;;
  esac
done

BASE=ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq
RELEASE=ftp://ftp.ncbi.nlm.nih.gov/refseq/release
SHARDS=$WORK/shards
LIST=$WORK/genomes.tsv                  # <url> <tab> <taxid>, one per assembly

if [ "$PLASMID" = 1 ] && [ "$SKIP_MAPS" = 1 ]; then
  echo "ERROR: plasmids have no assembly accession, so their taxids can only come" >&2
  echo "       from accession2taxid. Use --no-plasmid or drop --skip-maps." >&2
  exit 1
fi

for tool in kraken2-build dustmasker parallel wget; do
  command -v "$tool" >/dev/null || { echo "ERROR: $tool not on PATH" >&2; exit 1; }
done
[ -x "$SELF_DIR/mask_one.sh" ] || { echo "ERROR: $SELF_DIR/mask_one.sh missing" >&2; exit 1; }

mkdir -p "$DB" "$FASTAS" "$WORK" "$SHARDS"
echo "[k2db] db=$DB  k=$KMER l=$MINIMIZER s=$SPACES  max-db-size=$MAXDB  threads=$THREADS"
echo "[k2db] groups: $LIBS   viral=$VIRAL plasmid=$PLASMID"

# --- 1. taxonomy ------------------------------------------------------------
# The accession2taxid maps are only needed for sequences whose headers do not
# already carry a taxid - which here means the plasmid release files.
if [ ! -s "$DB/taxonomy/nodes.dmp" ]; then
  echo "[k2db] downloading taxonomy"
  if [ "$SKIP_MAPS" = 1 ]; then
    kraken2-build --db "$DB" --download-taxonomy --use-ftp --skip-maps
  else
    kraken2-build --db "$DB" --download-taxonomy --use-ftp
  fi
else
  echo "[k2db] taxonomy already present"
fi

# --- 2. UniVec_Core ---------------------------------------------------------
# Vectors, adapters and linkers. kraken2 tags these itself as taxid 81077
# (artificial sequences). They are a sink for lab artefacts: without them, an
# adapter read gets forced onto whatever real genome happens to share a stretch.
if [ ! -d "$DB/library/UniVec_Core" ]; then
  echo "[k2db] downloading UniVec_Core"
  kraken2-build --db "$DB" --download-library UniVec_Core --use-ftp --no-masking
fi

# --- 3. work out what to download -------------------------------------------
# assembly_summary columns: 5 = refseq_category, 6 = taxid, 20 = ftp_path.
# We keep only reference/representative assemblies, except for viruses, which
# are not flagged that way - there we take everything.
if [ ! -s "$LIST" ]; then
  : > "$LIST"
  for g in $LIBS; do
    echo "[k2db] assembly summary: $g"
    wget -q -O "$WORK/assembly_summary_$g.txt" "$BASE/$g/assembly_summary.txt"
    awk -F'\t' 'BEGIN{OFS="\t"}
      !/^#/ && ($5=="reference genome" || $5=="representative genome") {
        n=split($20,p,"/"); print $20"/"p[n]"_genomic.fna.gz", $6
      }' "$WORK/assembly_summary_$g.txt" >> "$LIST"
  done

  if [ "$VIRAL" = 1 ]; then
    echo "[k2db] assembly summary: viral (all assemblies)"
    wget -q -O "$WORK/assembly_summary_viral.txt" "$BASE/viral/assembly_summary.txt"
    awk -F'\t' 'BEGIN{OFS="\t"}
      !/^#/ { n=split($20,p,"/"); print $20"/"p[n]"_genomic.fna.gz", $6 }' \
      "$WORK/assembly_summary_viral.txt" >> "$LIST"
  fi

  [ "$TEST_N" -gt 0 ] && { head -n "$TEST_N" "$LIST" > "$LIST.tmp"; mv "$LIST.tmp" "$LIST"; }
fi
echo "[k2db] $(wc -l < "$LIST") assemblies to fetch"

# --- 4. download ------------------------------------------------------------
# -c resumes, so re-running the script does not re-fetch what is already there.
echo "[k2db] downloading genomes ($DL_JOBS parallel)"
cut -f1 "$LIST" | parallel -j "$DL_JOBS" wget -q -c -P "$FASTAS" {}

if [ "$PLASMID" = 1 ] && [ -z "$(ls "$FASTAS"/plasmid.*.genomic.fna.gz 2>/dev/null)" ]; then
  echo "[k2db] downloading plasmids"
  wget -q -c -r -nd -np -l1 -A 'plasmid.*.genomic.fna.gz' -P "$FASTAS" "$RELEASE/plasmid/"
fi

# --- 5. dustmask, tag with the taxid, add to the library --------------------
# One shard per parallel slot ({%}), so appends within a slot stay serial.
# add-to-library itself must run serially - kraken2 races on shared db state.
if [ -d "$DB/library/added" ] && [ -n "$(ls "$DB/library/added" 2>/dev/null)" ]; then
  echo "[k2db] library already populated, skipping masking"
else
  echo "[k2db] masking $(wc -l < "$LIST") genomes on $THREADS cores"
  awk -F'\t' -v d="$FASTAS" 'BEGIN{OFS="\t"} { n=split($1,p,"/"); print d"/"p[n], $2 }' "$LIST" \
    > "$WORK/tomask.tsv"
  parallel -j "$THREADS" --colsep '\t' \
    "$SELF_DIR/mask_one.sh {1} {2} $SHARDS/shard_{%}.fna" :::: "$WORK/tomask.tsv"

  # Plasmids carry no taxid we can inject, so pass "-" and let kraken2 resolve
  # them from accession2taxid during the build.
  if [ "$PLASMID" = 1 ]; then
    n_plasmid=$(find "$FASTAS" -maxdepth 1 -name 'plasmid.*.genomic.fna.gz' | wc -l)
    echo "[k2db] masking $n_plasmid plasmid files"
    [ "$n_plasmid" -gt 0 ] || echo "[k2db] WARN: no plasmid files found in $FASTAS" >&2
    find "$FASTAS" -maxdepth 1 -name 'plasmid.*.genomic.fna.gz' \
      | parallel -j "$THREADS" "$SELF_DIR/mask_one.sh {} - $SHARDS/plasmid_{%}.fna"
  fi

  for shard in "$SHARDS"/*.fna; do
    [ -s "$shard" ] || continue
    echo "[k2db]   add-to-library $(basename "$shard") ($(du -h "$shard" | cut -f1))"
    kraken2-build --db "$DB" --add-to-library "$shard" --no-masking
    rm -f "$shard"
  done
fi

# --- 6. build ---------------------------------------------------------------
# --minimizer-spaces must be given explicitly: kraken2-build defaults to 7, and
# every masked position costs a base of specificity (matching uses l - s bases).
# --no-masking because we dustmasked already.
echo "[k2db] building index"
kraken2-build --build --db "$DB" \
  --kmer-len "$KMER" \
  --minimizer-len "$MINIMIZER" \
  --minimizer-spaces "$SPACES" \
  --max-db-size "$MAXDB" \
  --threads "$THREADS" \
  --no-masking

echo "[k2db] done: $DB"
ls -lh "$DB"/*.k2d
echo "[k2db] check the header before using it:"
echo "       kraken2-inspect --db $DB | head -7"

if [ -n "$STORAGE" ]; then
  echo "[k2db] copying runtime database to $STORAGE"
  mkdir -p "$STORAGE"
  rsync -a --exclude 'library/' "$DB"/ "$STORAGE"/
fi

