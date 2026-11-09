#!/bin/bash
# AF3 data pipeline ONLY (MSA search; template search is off via "templates": []).
# CPU node, no GPU. Writes msas/af3/<name>/<name>_data.json, which is reused for every seed.
#   sbatch scripts/af3_msa.slurm <TARGET>
#SBATCH --job-name=af3msa
#SBATCH --partition=general
#SBATCH --account=dkhf3-lab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=/cluster/pixstor/dkhf3-lab/users/nancy/PROJECT4/logs/af3msa_%x_%j.out

set -o pipefail
shopt -s expand_aliases
T=${1:?target}
P=/cluster/pixstor/dkhf3-lab/users/nancy/PROJECT4
OUT=$P/msas/af3
mkdir -p "$OUT"
LOG=$OUT/${T}_msa.log
export TMPDIR=$OUT/tmp_${SLURM_JOB_ID}; mkdir -p "$TMPDIR"
module load alphafold3/alphafold3_v301_deepmind

echo "[$(date)] $T on $(hostname), $(nproc) CPUs"
T0=$(date +%s.%N)
run_script_a3 --json_path=$P/data/inputs/af3_msa/$T.json \
    --db_dir=/cluster/VAST/cietest/alphafold3_database \
    --model_dir=/cluster/VAST/cietest/alphafold3_model_parameters \
    --output_dir="$OUT" --norun_inference \
    --jackhmmer_n_cpu=8 --nhmmer_n_cpu=8 2>&1 | tee "$LOG"
rc=$?
wall=$(echo "$(date +%s.%N) - $T0" | bc)
pipe=$(grep -oE "Running data pipeline for chain .* took [0-9.]+ seconds" "$LOG" | grep -oE "[0-9.]+ seconds" | awk '{s+=$1} END{print s+0}')
printf "%s\taf3\t%s\t%s\t%.1f\t%s\n" "$T" "$SLURM_JOB_ID" "$rc" "$wall" "$pipe" >> $P/msas/msa_timing_af3.tsv
rm -rf "$TMPDIR"
exit $rc
