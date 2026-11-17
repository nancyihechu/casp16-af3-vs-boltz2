#!/bin/bash
# DIAGNOSTIC Boltz-2 run (never part of the benchmark tables).
#   sbatch scripts/run_boltz2_diag.slurm <LABEL> <TARGET> <SEED> [extra boltz flags...]
# Same settings as scripts/run_boltz2.slurm (cached MSA, no templates, SMILES only,
# --diffusion_samples 5), plus any extra flags. Output: predictions/diagnostics/boltz2_<LABEL>/<T>/s<SEED>
#
# Potentials in boltz 2.2.1: the only CLI flag is --use_potentials ("Whether to use potentials
# for steering. Default is False."), which sets fk_steering and physical_guidance_update.
# Potentials are therefore OFF unless that flag is passed; the benchmark runs never passed it.
#SBATCH --job-name=bz_diag
#SBATCH --partition=dkhf3-lab-gpu
#SBATCH --account=dkhf3-lab
#SBATCH --gres=gpu:L40S:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/cluster/pixstor/dkhf3-lab/users/nancy/PROJECT4/logs/boltz2_diag_%j.out

set -o pipefail
L=${1:?label}; T=${2:?target}; S=${3:?seed}
P=/cluster/pixstor/dkhf3-lab/users/nancy/PROJECT4
OUT=$P/predictions/diagnostics/boltz2_$L/$T/s$S
mkdir -p "$OUT"
export BOLTZ_CACHE=/cluster/pixstor/dkhf3-lab/users/nancy/envs/boltz2_cache
BOLTZ=/cluster/pixstor/dkhf3-lab/users/nancy/envs/boltz2/bin/boltz

echo "[$(date)] DIAGNOSTIC $L $T seed $S extra: ${*:4}  on $(hostname)"
T0=$(date +%s.%N)
$BOLTZ predict "$P/data/inputs/boltz/$T.yaml" --out_dir "$OUT" --cache "$BOLTZ_CACHE" \
    --model boltz2 --seed "$S" --diffusion_samples 5 --output_format mmcif --override "${@:4}" 2>&1 | tee "$OUT/boltz.log"
rc=$?
wall=$(echo "$(date +%s.%N) - $T0" | bc)
printf "label\ttarget\tseed\tjob_id\trc\twall_s\textra_flags\n%s\t%s\t%s\t%s\t%s\t%.1f\t%s\n" \
    "$L" "$T" "$S" "$SLURM_JOB_ID" "$rc" "$wall" "${*:4}" > "$OUT/timing.tsv"
exit $rc
