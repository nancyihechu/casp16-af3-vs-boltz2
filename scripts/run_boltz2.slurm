#!/bin/bash
# Boltz-2 (envs/boltz2, boltz 2.2.1) for one target and one seed.
#   sbatch scripts/run_boltz2.slurm <TARGET> <SEED>
# Input: data/inputs/boltz/<TARGET>.yaml with cached MSA paths (msa: field).
# NO --use_msa_server: if a cached MSA is missing, Boltz errors out instead of silently
# fetching a new one. Templates: none given in the YAML, so off. Ligands are SMILES only.
# --diffusion_samples 5 per seed, matching AF3's 5 samples per seed.
# All other settings are Boltz defaults (recycling_steps 3, sampling_steps 200, no potentials).
#SBATCH --job-name=boltz2
#SBATCH --partition=dkhf3-lab-gpu
#SBATCH --account=dkhf3-lab
#SBATCH --gres=gpu:L40S:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/cluster/pixstor/dkhf3-lab/users/nancy/PROJECT4/logs/boltz2_%j.out

set -o pipefail
T=${1:?target}; S=${2:?seed}
P=/cluster/pixstor/dkhf3-lab/users/nancy/PROJECT4
OUT=$P/predictions/boltz2/$T/s$S
mkdir -p "$OUT"
export BOLTZ_CACHE=/cluster/pixstor/dkhf3-lab/users/nancy/envs/boltz2_cache
BOLTZ=/cluster/pixstor/dkhf3-lab/users/nancy/envs/boltz2/bin/boltz

echo "[$(date)] $T seed $S on $(hostname) GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
T0=$(date +%s.%N)
$BOLTZ predict "$P/data/inputs/boltz/$T.yaml" --out_dir "$OUT" --cache "$BOLTZ_CACHE" \
    --model boltz2 --seed "$S" --diffusion_samples 5 --output_format mmcif --override 2>&1 | tee "$OUT/boltz.log"
rc=$?
wall=$(echo "$(date +%s.%N) - $T0" | bc)
printf "name\tjob_id\trc\twall_s\tdata_pipeline_s\tinference_s\tmsa_source\n%s_s%s\t%s\t%s\t%.1f\t0\t%.1f\t%s\n" \
    "$T" "$S" "$SLURM_JOB_ID" "$rc" "$wall" "$wall" "ColabFold MMseqs2 server via boltz.main.compute_msa (cached csv)" > "$OUT/timing.tsv"
cat "$OUT/timing.tsv"
exit $rc
