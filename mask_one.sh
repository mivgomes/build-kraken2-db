#!/usr/bin/env bash
# Dustmask one genome and append it to a shard. Called by GNU parallel.
#
# Usage: mask_one.sh <genome.fna.gz> <taxid|-> <output_shard>
#
#   dustmasker soft-masks low-complexity regions; awk then hard-masks them
#   (lowercase -> N) and, unless the taxid is "-", injects >kraken:taxid|<taxid>|
#   into every header. Pass "-" for sequences that have no assembly-level taxid
#   (the plasmid release files); kraken2 resolves those from accession2taxid.
#
#   A genome that fails is logged and skipped (exit 0) so one bad file does not
#   kill the whole parallel run.

set -uo pipefail

f="$1"; t="$2"; shard="$3"

if [ "$t" = "-" ]; then
  tag='1'                                    # keep headers as they are
else
  tag='/^>/{sub(/^>/,">kraken:taxid|"t"|")} 1'
fi

if ! { zcat "$f" \
        | "${DUSTMASKER:-dustmasker}" -infmt fasta -outfmt fasta \
        | awk -v t="$t" '!/^>/{gsub(/[a-z]/,"N")} '"$tag" \
        >> "$shard"; }; then
  echo "[k2db] WARN: masking failed for $(basename "$f") -- skipping" >&2
fi
exit 0
