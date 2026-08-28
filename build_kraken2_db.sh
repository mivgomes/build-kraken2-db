#!/usr/bin/env bash
# Build a custom Kraken2 database (k=22) for ancient-DNA metagenomic screening.
#
# WHY: the prebuilt long-k DBs (k2_core_nt k=35, etc.) leave ~half of hg19-confirmed
# human reads unclassified, because aDNA reads are short and C->T damaged. k=22
# recovers them. This is the Kraken2 twin of the KrakenUniq k=22 build, so the
# benchmark can separate the k-mer effect from the tool effect. Built from the same
# RefSeq release-232 genomes (Vernot et al. 2021 sediment-DNA recipe).
#
# This runs by a SLURM wrapper (build_kraken2_db.sbatch). Everything is env-driven;
# Progressive deletion saves space.
#
# Same inputs / masking / sharding as the KrakenUniq build; the differences are:
#   - kraken2-build instead of krakenuniq-build
#   - the index build uses minimizers (no jellyfish); --no-masking because we pre-mask
#   - none of the KrakenUniq workarounds (no libexec PATH hack, no exit-255 handling)

set -euo pipefail

# Self-diagnosing error trap: if any command fails, print the exact failing line +
# context instead of dying silently.
err_report() {
  local exit_code=$?
  echo "[k2db] *** FAILURE on line ${1} (exit ${exit_code}) ***" >&2
  echo "[k2db] command: ${BASH_COMMAND}" >&2
  echo "[k2db] PATH=${PATH}" >&2
  exit "$exit_code"
}
trap 'err_report ${LINENO}' ERR

# paths (inputs already copied to scratch)
# This script is meant to live IN the build folder (e.g. /SCRATCH/build-kraken2-db/)
# next to fastas/ and assembly_summary_*.txt, so SCR defaults to the script's own
# directory. Override SCR to point elsewhere.
SELF_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCR=${SCR:-$SELF_DIR} # path to where RefSeq_Release232/* was copied
FASTAS_DIR=${FASTAS_DIR:-$SCR/fastas} # *_genomic.fna.gz (24,031 files, no viral)
SUMMARY_DIR=${SUMMARY_DIR:-$SCR} # assembly_summary_<group>.txt
DB=${DB:-$SCR/kraken2_refseq232_k22} # DB built here (fast local disk)
SHARD_DIR=${SHARD_DIR:-$DB.shards} # transient masked shards (deleted as added)

# build parameters
# Kraken2 minimizer constraints: MINIMIZER (l) must be <= KMER (k) and <= 31;
# SPACES (s) must be small relative to l (kraken2-build validates and errors with the
# exact max). For k=22 we use l=20, s=5. Spaced seeds also tolerate mismatches, which
# helps C->T-damaged reads.
KMER=${KMER:-22}
MINIMIZER=${MINIMIZER:-20}
SPACES=${SPACES:-5}
THREADS=${THREADS:-8}
NSHARDS=${NSHARDS:-16} # masked genomes are concatenated into this many shards

# run modes
PILOT_N=${PILOT_N:-0} # >0 : only the first N genomes (de-risk time/size + taxid headers)
DELETE_GZ=${DELETE_GZ:-0} # 1  : delete the staged *.gz after masking (original stays on /data)
STORAGE_DIR=${STORAGE_DIR:-} # set: rsync the finished runtime DB here (e.g. /eq_Peyregne/.../kraken2_refseq232_k22)

# tools (override if not on PATH)
K2BUILD=${K2BUILD:-kraken2-build}
DUSTMASKER=${DUSTMASKER:-dustmasker}

echo "[k2db] SCR=$SCR"
echo "[k2db] DB=$DB  k=$KMER  l=$MINIMIZER  s=$SPACES  threads=$THREADS  shards=$NSHARDS"
[ "$PILOT_N" -gt 0 ] && echo "[k2db] *** PILOT MODE: first $PILOT_N genomes ***"

# 0) sanity: inputs present
n_fastas=$(printf '%s\n' "$FASTAS_DIR"/*_genomic.fna.gz 2>/dev/null | wc -l)
n_sum=$(printf '%s\n' "$SUMMARY_DIR"/assembly_summary_*.txt 2>/dev/null | wc -l)
if [ "$n_fastas" -eq 0 ] || [ "$n_sum" -eq 0 ]; then
  echo "ERROR: need *_genomic.fna.gz in $FASTAS_DIR and assembly_summary_*.txt in $SUMMARY_DIR" >&2
  echo "       found $n_fastas genomes / $n_sum summaries. Set SCR to your staged copy." >&2
  exit 1
fi
command -v "$K2BUILD"    >/dev/null || { echo "ERROR: $K2BUILD not on PATH (conda activate taxclassification)" >&2; exit 1; }
command -v "$DUSTMASKER" >/dev/null || { echo "ERROR: $DUSTMASKER not on PATH" >&2; exit 1; }
echo "[k2db] inputs OK: $n_fastas genomes, $n_sum summaries"

mkdir -p "$DB" "$SHARD_DIR"

# 1) Extract taxonomy files (names.dmp / nodes.dmp) from NCBI
# (--skip-maps could be added: taxids come from injected headers, so the multi-GB
#  accession2taxid maps are not strictly needed; plain download is the safe default.)
if [ ! -s "$DB/taxonomy/nodes.dmp" ]; then
  echo "[k2db] downloading taxonomy"
  "$K2BUILD" --db "$DB" --download-taxonomy --use-ftp --skip-maps
else
  echo "[k2db] taxonomy already present, skipping download"
fi

# 2) accession -> taxid map (8 non-viral summaries; col1=acc col6=taxid)
declare -A TAXID
while IFS=$'\t' read -r acc t; do
  [ -n "$acc" ] && TAXID["$acc"]="$t"
done < <(awk -F '\t' 'FNR>1 && $0 !~ /^#/ {print $1"\t"$6}' \
             "$SUMMARY_DIR"/assembly_summary_*.txt)
echo "[k2db] taxid map: ${#TAXID[@]} accessions"


# 3) mask + hard-mask + tag + shard + add-to-library
# Per genome:
#   dustmasker identifies repetitive/low-complexity regions and we hard-mask them to N
#   Adds >kraken:taxid|<taxid>| into every header, so the sequences are tagged
#   Genomes are concatenated into NSHARDS shards; each shard is added to the library
#   then deleted, so disk never holds the whole masked FASTA at once.
# Identical to the KrakenUniq masking -- keeping it the same makes the two DBs comparable.

# If the library is already populated, skip masking/add so the index build (step 4)
# can be re-run without redoing the expensive masking step.
if [ -d "$DB/library/added" ] && \
   [ "$(find "$DB/library/added" -type f | wc -l)" -gt 0 ]; then
    echo "[k2db] library already populated; skipping masking/add-to-library"
    SKIP_LIBRARY=1
else
    SKIP_LIBRARY=0
fi

if [ "$SKIP_LIBRARY" -eq 0 ]; then
  mapfile -t FASTAS < <(printf '%s\n' "$FASTAS_DIR"/*_genomic.fna.gz | sort)
  [ "$PILOT_N" -gt 0 ] && FASTAS=("${FASTAS[@]:0:$PILOT_N}")
  total=${#FASTAS[@]}
  per_shard=$(( (total + NSHARDS - 1) / NSHARDS )); [ "$per_shard" -lt 1 ] && per_shard=1

  # adds current shard to the library, then frees its space
  flush_shard() {
    local sf="$1" # shard file
    [ -s "$sf" ] || { rm -f "$sf"; return; }
    echo "[k2db] add-to-library $(basename "$sf") ($(du -h "$sf" | cut -f1))"
    # kraken2-build --add-to-library exits non-zero on real failure (the ERR trap
    # catches it) -- no exit-255 workaround needed, unlike KrakenUniq. Note Kraken2
    # stores the file under library/added/ with a generated name, so we don't check
    # for the original basename; we trust the exit status.
    "$K2BUILD" --db "$DB" --add-to-library "$sf" --no-masking
    rm -f "$sf" # free space
  }

  shard=0; in_shard=0; masked=0; missing=0
  shard_file="$SHARD_DIR/shard_$shard.fna"; : > "$shard_file"
  for f in "${FASTAS[@]}"; do
    bn=$(basename "$f")
    acc=$(printf '%s' "$bn" | grep -oE '^GC[FA]_[0-9]+\.[0-9]+' || true) # accession ID
    t="${TAXID[$acc]:-}" # taxid for this accession (from the summary files)
    if [ -z "$t" ]; then
      echo "[k2db] WARN: no taxid for '$acc' ($bn) -- skipping" >&2
      missing=$((missing+1)); continue
    fi

    # mask + inject >kraken:taxid|<taxid>| header + hard-mask (lowercase -> N)
    if ! { zcat "$f" \
            | "$DUSTMASKER" -infmt fasta -outfmt fasta \
            | awk -v t="$t" '/^>/{sub(/^>/,">kraken:taxid|"t"|")} !/^>/{gsub(/[a-z]/,"N")} 1' \
            >> "$shard_file"; }; then
      echo "[k2db] WARN: masking failed for $bn -- skipping" >&2
      continue
    fi

    masked=$((masked+1)); in_shard=$((in_shard+1))
    if (( masked % 500 == 0 )); then echo "[k2db] masked $masked / $total"; fi

    if [ "$in_shard" -ge "$per_shard" ]; then
      flush_shard "$shard_file"
      shard=$((shard+1)); in_shard=0
      shard_file="$SHARD_DIR/shard_$shard.fna"; : > "$shard_file"
    fi
  done

  flush_shard "$shard_file" # last partial shard
  rmdir "$SHARD_DIR" 2>/dev/null || true
  echo "[k2db] masked $masked genomes into library ($missing skipped for missing taxid)"

  # Optional: clean temporary .gz (the original should live on /path/to/storage/folder/).
  if [ "$DELETE_GZ" = 1 ] && [ "$PILOT_N" -eq 0 ]; then
    echo "[k2db] deleting staged gz to free space"
    rm -f "$FASTAS_DIR"/*_genomic.fna.gz
  fi

fi

# 4) build the index (minimizer-based; lighter/faster than KrakenUniq's exact k-mers)
echo "[k2db] building index (k=$KMER, l=$MINIMIZER, s=$SPACES, threads=$THREADS)"
# --no-masking: we already hard-masked above, so don't let kraken2-build dustmask again.
"$K2BUILD" --db "$DB" --build \
  --kmer-len "$KMER" --minimizer-len "$MINIMIZER" --minimizer-spaces "$SPACES" \
  --threads "$THREADS" --no-masking

echo "[k2db] BUILD DONE: $DB"
ls -lh "$DB"/*.k2d 2>/dev/null || true
echo "[k2db] verify with:  kraken2-inspect --db $DB | grep -E 'Homo sapiens|Mammalia'"

# Optional: move DB to /path/to/storage/folder/. Only the three *.k2d files are needed
# at classification time (kraken2-build --clean removes library/ + taxonomy/ if you
# want to slim the local copy first).
if [ -n "$STORAGE_DIR" ]; then
  echo "[k2db] persisting runtime DB -> $STORAGE_DIR"
  mkdir -p "$STORAGE_DIR"
  rsync -a --exclude 'library/' "$DB"/ "$STORAGE_DIR"/
  echo "[k2db] persisted. Point config.yaml taxonomy.databases.kraken2 at: $STORAGE_DIR"
fi
