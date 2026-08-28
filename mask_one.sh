#!/usr/bin/env bash
# Mask + hard-mask + taxid-tag ONE genome, for the parallel step of
# build_kraken2_db.sh. Runs as a standalone script (called by GNU parallel) rather
# than an exported bash function, which is the robust way to use parallel.
#
# Usage: mask_one.sh <genome.fna.gz> <taxid> <output_shard>   (appends to the shard)
#   - dustmasker soft-masks low-complexity -> awk hard-masks (lowercase -> N) and
#     injects >kraken:taxid|<taxid>| into every header.
#   - A single genome that fails is logged and skipped (exit 0), so it never aborts
#     the whole parallel run.
# DUSTMASKER may be set in the environment (inherited from the parent); defaults to
# `dustmasker` on PATH.
set -uo pipefail

f="$1"; t="$2"; shard="$3"

if ! { zcat "$f" \
        | "${DUSTMASKER:-dustmasker}" -infmt fasta -outfmt fasta \
        | awk -v t="$t" '/^>/{sub(/^>/,">kraken:taxid|"t"|")} !/^>/{gsub(/[a-z]/,"N")} 1' \
        >> "$shard"; }; then
  echo "[k2db] WARN: masking failed for $(basename "$f") -- skipping" >&2
fi
exit 0
