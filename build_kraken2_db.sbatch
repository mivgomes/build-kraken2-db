#!/usr/bin/env bash
#SBATCH --job-name=k2db_build
#SBATCH --account=eq_peyregne
#SBATCH --partition=debug
#SBATCH --nodelist=geddy
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --hint=compute_bound
#SBATCH --output=/data/home/mgomes/logs/k2db_build_%j.out
#SBATCH --error=/data/home/mgomes/logs/k2db_build_%j.err
#
# Create logs folder on local disk first:  mkdir -p ~/logs
# Both files live IN the build folder (e.g. /scratch2/build-kraken2-db/) next to
# fastas/ + assembly_summary_*.txt.
#
#   cd /scratch2/mgomes/build-kraken2-db
#
# Test run first (~50 genomes incl. human GRCh38 + mouse, throwaway DB):
#   PILOT_N=50 DB=$PWD/DB_pilot sbatch build_kraken2_db.sbatch
#
# Full run (after the pilot is OK), optionally copying the result to storage:
#   DELETE_GZ=1 STORAGE_DIR=/data/home/mgomes/data/ines/kraken2_refseq232_k22 sbatch build_kraken2_db.sbatch

# NOTE: the masking (step 3) is the serial dustmasker loop as the
# KrakenUniq build, so it is still the bottleneck (~days for the full set). Only the
# index build (step 4) is lighter/faster for Kraken2.

set -euo pipefail

# conda env with kraken2 / kraken2-build (self-contained -- no libexec PATH hack, no
# jellyfish, unlike KrakenUniq).
export PATH=/data/home/mgomes/.conda/envs/blast_env/bin:$PATH
export THREADS=${SLURM_CPUS_PER_TASK:-8}
export OMP_NUM_THREADS=${THREADS}

# The build folder = where we submitted from (SLURM_SUBMIT_DIR), or override BUILD_DIR.
BUILD_DIR=${BUILD_DIR:-${SLURM_SUBMIT_DIR:-/scratch2/mgomes/build-kraken2-db}}
cd "$BUILD_DIR"
echo "[sbatch] host = $(hostname) job = ${SLURM_JOB_ID:-none} threads = $THREADS"
echo "[sbatch] build dir: $BUILD_DIR"
echo "[sbatch] kraken2-build: $(command -v kraken2-build)"
echo "[sbatch] dustmasker: $(command -v dustmasker)"
# --mem=700G sets SLURM_MEM_PER_NODE (MB); SLURM_MEM_PER_CPU stays unset in that case.
echo "[sbatch] requested cpus=$THREADS  mem=${SLURM_MEM_PER_NODE:-$(( ${SLURM_MEM_PER_CPU:-0} * ${SLURM_CPUS_PER_TASK:-1} ))}MB"

# Run the build script
bash "$BUILD_DIR/build_kraken2_db.sh"
