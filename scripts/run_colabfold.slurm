#!/bin/bash
# AF2-Multimer v3 via ColabFold 1.5.5 (envs/colabfold) for one target.
#   sbatch scripts/run_colabfold.slurm <TARGET> <FIRST_SEED> <NUM_SEEDS>
# Seeds run are range(FIRST_SEED, FIRST_SEED+NUM_SEEDS); each seed runs all 5 models,
# so there are 5 predictions per seed.
# Input: the cached msas/colabfold/<T>.a3m, so the MSA server is not called again.
# Templates: no --templates flag, so off. No Amber relax. Other settings are ColabFold defaults.
#SBATCH --job-name=af2m
#SBATCH --partition=dkhf3-lab-gpu
#SBATCH --account=dkhf3-lab
#SBATCH --gres=gpu:L40S:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/cluster/pixstor/dkhf3-lab/users/nancy/PROJECT4/logs/af2m_%j.out

set -o pipefail
T=${1:?target}; S=${2:?first seed}; N=${3:?num seeds}
P=/cluster/pixstor/dkhf3-lab/users/nancy/PROJECT4
OUT=$P/predictions/af2m/$T/s${S}_n${N}
mkdir -p "$OUT"

echo "[$(date)] $T seeds $S..$((S+N-1)) on $(hostname) GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
T0=$(date +%s.%N)
/cluster/pixstor/dkhf3-lab/users/nancy/envs/colabfold/bin/colabfold_batch \
    "$P/msas/colabfold/$T.a3m" "$OUT" \
    --model-type alphafold2_multimer_v3 --num-models 5 \
    --random-seed "$S" --num-seeds "$N" \
    --data /cluster/pixstor/dkhf3-lab/users/nancy/colabfold_weights 2>&1 | tee "$OUT/colabfold.log"
rc=$?
wall=$(echo "$(date +%s.%N) - $T0" | bc)
printf "name\tjob_id\trc\twall_s\tdata_pipeline_s\tinference_s\tmsa_source\n%s_s%s_n%s\t%s\t%s\t%.1f\t0\t%.1f\t%s\n" \
    "$T" "$S" "$N" "$SLURM_JOB_ID" "$rc" "$wall" "$wall" "ColabFold MMseqs2 server (cached a3m)" > "$OUT/timing.tsv"
cat "$OUT/timing.tsv"
exit $rc
