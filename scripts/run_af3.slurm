#!/bin/bash
# AlphaFold3 (DeepMind v3.0.1 container, git v3.0.1-13-ga3cf058) for one input JSON.
#   sbatch scripts/run_af3.slurm <input.json> <output_dir> [extra run_alphafold flags...]
#   Per-seed runs pass the cached <name>_data.json and --norun_data_pipeline.
# Templates must already be disabled in the JSON ("templates": [] on every protein
# chain, no MSA fields). MSA comes from AF3's own jackhmmer/nhmmer pipeline.
# No automatic OOM retry: a failure is reported, not worked around.
#SBATCH --job-name=af3
#SBATCH --partition=dkhf3-lab-gpu
#SBATCH --account=dkhf3-lab
#SBATCH --gres=gpu:L40S:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/cluster/pixstor/dkhf3-lab/users/nancy/PROJECT4/logs/af3_%j.out

set -o pipefail
shopt -s expand_aliases

INPUT=$(readlink -f "${1:?usage: sbatch run_af3.slurm <input.json> <output_dir>}")
OUT=${2:?usage: sbatch run_af3.slurm <input.json> <output_dir>}
mkdir -p "$OUT"; OUT=$(readlink -f "$OUT")
NAME=$(basename "$INPUT" .json)
LOG=$OUT/${NAME}_af3.log
TIMING=$OUT/${NAME}_timing.tsv

export TMPDIR=$OUT/tmp_${SLURM_JOB_ID}
mkdir -p "$TMPDIR"

module load alphafold3/alphafold3_v301_deepmind

echo "[$(date)] Extra flags: ${*:3}"
echo "[$(date)] Input: $INPUT  Host: $(hostname)  GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader)"
T0=$(date +%s.%N)
run_script_a3 \
    --json_path="$INPUT" \
    --model_dir=/cluster/VAST/cietest/alphafold3_model_parameters \
    --db_dir=/cluster/VAST/cietest/alphafold3_database \
    --output_dir="$OUT" \
    --jackhmmer_n_cpu=8 --nhmmer_n_cpu=8 "${@:3}" 2>&1 | tee "$LOG"
rc=$?
T1=$(date +%s.%N)

# Split wall-clock into MSA/data pipeline vs inference, from AF3's own log lines.
data_s=$(grep -oE "Running data pipeline for chain .* took [0-9.]+ seconds" "$LOG" | grep -oE "[0-9.]+ seconds" | awk '{s+=$1} END{print s+0}')
infer_s=$(grep -oE "Running model inference with seed [0-9]+ took [0-9.]+ seconds" "$LOG" | grep -oE "took [0-9.]+" | awk '{s+=$2} END{print s+0}')
wall_s=$(echo "$T1 - $T0" | bc)
printf "name\tjob_id\trc\twall_s\tdata_pipeline_s\tinference_s\tmsa_source\n%s\t%s\t%s\t%.1f\t%s\t%s\t%s\n" \
    "$NAME" "$SLURM_JOB_ID" "$rc" "$wall_s" "$data_s" "$infer_s" "AF3 jackhmmer/nhmmer local (/cluster/VAST/cietest/alphafold3_database)" > "$TIMING"
cat "$TIMING"

echo "[$(date)] AlphaFold3 returned $rc"
rm -rf "$TMPDIR"
exit $rc
