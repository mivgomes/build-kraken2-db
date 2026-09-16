#!/usr/bin/env bash
#SBATCH --job-name=k2db
#SBATCH --account=eq_peyregne
#SBATCH --partition=debug
#SBATCH --nodelist=geddy
#SBATCH --cpus-per-task=32
#SBATCH --mem=1200G
#SBATCH --hint=compute_bound
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#
# Everything after the script name is passed straight to build_db.sh, so the
# parameters live on the command line:
#
#   sbatch build_db.sbatch --db /scratch-peyregne/mgomes/k2_refseq_k25 \
#                          --fastas /scratch-peyregne/mgomes/fastas \
#                          --storage /data/home/mgomes/data/ines/k2_refseq_k25
#
# Try the recipe on a handful of genomes first:
#
#   sbatch build_db.sbatch --db $PWD/db_test --test 50 --no-plasmid --skip-maps
#
# --mem must stay above --max-db-size (1.1 TB = 1024 GiB) with room to spare.

set -euo pipefail

export PATH=/data/home/mgomes/.conda/envs/blast_env/bin:$PATH
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

cd "${SLURM_SUBMIT_DIR:-$PWD}"

echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID:-none} cpus=${SLURM_CPUS_PER_TASK:-?}"
echo "[sbatch] kraken2-build: $(command -v kraken2-build)"
echo "[sbatch] dustmasker:    $(command -v dustmasker)"
echo "[sbatch] parallel:      $(command -v parallel)"

bash build_db.sh --threads "${SLURM_CPUS_PER_TASK:-8}" "$@"
