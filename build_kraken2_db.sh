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
KMER=${KMER:-22}
MINIMIZER=${MINIMIZER:-22}
THREADS=${THREADS}          # parallel masking jobs (one dustmasker per core) + build threads
ROUND_SIZE=${ROUND_SIZE:-2000} # genomes masked per round before add-to-library (bounds scratch peak)

# run modes
TEST_N=${TEST_N:-0} # >0 : only the first N genomes (de-risk time/size + taxid headers)
DELETE_GZ=${DELETE_GZ:-0} # 1  : delete the staged *.gz after masking (original stays on /data)
STORAGE_DIR=${STORAGE_DIR:-} # set: rsync the finished runtime DB here (e.g. /eq_Peyregne/.../kraken2_refseq232_k22)

# tools (override if not on PATH)
K2BUILD=${K2BUILD:-kraken2-build}
DUSTMASKER=${DUSTMASKER:-dustmasker}
PARALLEL=${PARALLEL:-parallel}   # GNU parallel (from your conda env; must be on PATH)

echo "[k2db] SCR=$SCR"
echo "[k2db] DB=$DB  k=$KMER  l=$MINIMIZER  s=0  threads=$THREADS  round=$ROUND_SIZE"
[ "$TEST_N" -gt 0 ] && echo "[k2db] *** TEST MODE: first $TEST_N genomes ***"

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
command -v "$PARALLEL"   >/dev/null || { echo "ERROR: $PARALLEL (GNU parallel) not on PATH -- activate the conda env that has it" >&2; exit 1; }
[ -x "$SELF_DIR/mask_one.sh" ] || { echo "ERROR: helper $SELF_DIR/mask_one.sh missing or not executable" >&2; exit 1; }
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


# 3) mask + hard-mask + tag + add-to-library  (PARALLEL masking)
# Per genome: dustmasker finds low-complexity regions -> mask_one.sh hard-masks them to
# N and injects >kraken:taxid|<taxid>| into every header. Masking is the slow, one-core-
# per-genome step, so we run it with GNU parallel across THREADS cores. Work is done in
# ROUND_SIZE-genome rounds: each round masks in parallel into per-slot shards ({%} = slot
# 1..THREADS, so appends within a slot are serial and never collide), then adds those
# shards to the library and deletes them -- bounding the scratch peak. add-to-library
# itself stays SERIAL (kraken2 races on shared DB state if run concurrently). Same masking
# recipe as KrakenUniq, so the two DBs stay comparable.

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
  [ "$TEST_N" -gt 0 ] && FASTAS=("${FASTAS[@]:0:$TEST_N}")

  # Phase A (fast, serial): build the work list  fasta<TAB>taxid, dropping any genome
  # whose accession has no taxid in the summaries.
  work="$SHARD_DIR/work.tsv"; : > "$work"
  missing=0
  for f in "${FASTAS[@]}"; do
    acc=$(printf '%s' "$(basename "$f")" | grep -oE '^GC[FA]_[0-9]+\.[0-9]+' || true)
    t="${TAXID[$acc]:-}"
    if [ -z "$t" ]; then
      echo "[k2db] WARN: no taxid for '$acc' -- skipping" >&2
      missing=$((missing+1)); continue
    fi
    printf '%s\t%s\n' "$f" "$t" >> "$work"
  done
  n_work=$(wc -l < "$work")
  echo "[k2db] to mask: $n_work genomes ($missing skipped for missing taxid); ${THREADS}-way parallel, round=$ROUND_SIZE"

  if [ "$n_work" -gt 0 ]; then
    export DUSTMASKER   # mask_one.sh reads it (env inherited by parallel's children)
    round=0
    split -l "$ROUND_SIZE" "$work" "$SHARD_DIR/round_"
    for chunk in "$SHARD_DIR"/round_*; do
      round=$((round+1))
      echo "[k2db] round $round: masking $(wc -l < "$chunk") genomes"
      # {1}=fasta {2}=taxid {%}=slot(1..THREADS) -> one shard per slot (no write clashes)
      "$PARALLEL" --jobs "$THREADS" --colsep '\t' \
        "$SELF_DIR/mask_one.sh {1} {2} $SHARD_DIR/shard_{%}.fna" :::: "$chunk"
      # add-to-library the round's shards SERIALLY, then free their space
      for shard in "$SHARD_DIR"/shard_*.fna; do
        [ -s "$shard" ] || { rm -f "$shard" 2>/dev/null || true; continue; }
        echo "[k2db]   add-to-library $(basename "$shard") ($(du -h "$shard" | cut -f1))"
        "$K2BUILD" --db "$DB" --add-to-library "$shard" --no-masking
        rm -f "$shard"
      done
      rm -f "$chunk"
    done
  fi

  rm -f "$work"
  rmdir "$SHARD_DIR" 2>/dev/null || true
  echo "[k2db] masking + add-to-library complete: $n_work genomes"

  # Optional: clean temporary .gz (the original should live on /path/to/storage/folder/).
  if [ "$DELETE_GZ" = 1 ] && [ "$TEST_N" -eq 0 ]; then
    echo "[k2db] deleting staged gz to free space"
    rm -f "$FASTAS_DIR"/*_genomic.fna.gz
  fi

fi

# 4) build the index (minimizer-based; lighter/faster than KrakenUniq's exact k-mers)
echo "[k2db] building index (k=$KMER, l=$MINIMIZER, s=0, threads=$THREADS)"
# --minimizer-spaces MUST be passed explicitly: kraken2-build defaults to 7, NOT 0.
# With l=22, s=7 exceeds the s <= l/4 ceiling (kraken2 rejects it), and any s>0
# shrinks the match to l-s informative bases -- the saturated-database failure mode.
# --no-masking: we already hard-masked above, so don't let kraken2-build dustmask again.
"$K2BUILD" --build --db "$DB" --kmer-len "$KMER" --minimizer-len "$MINIMIZER" --minimizer-spaces 0 --threads "$THREADS" --no-masking --max-db-size 1000000000000

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
