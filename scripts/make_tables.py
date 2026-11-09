"""Build results/tables/timing.csv and results/tables/confidence.csv from the run manifest.

  envs/p3_tools/bin/python -I scripts/make_tables.py   (run scripts/collect_runs.py first)

timing.csv: one row per (target, method, seed) run.
  msa_s and inference_s are SEPARATE columns and are never summed.
  msa_s is the one-time per-target MSA pass (msa_s_scope = per_target_once). It is
  repeated on each seed row for convenience; it was NOT paid once per seed.
  msa_source:
    af3_local_jackhmmer       AF3 jackhmmer/nhmmer on local DBs (CPU node, a real search)
    colabfold_server_cached   ColabFold MMseqs2 server returned a cached result
    colabfold_server_fresh    server ran a new search
    colabfold_server_mixed    some requests cached, some fresh (see msas/msa_server_cache_status.tsv)
    none                      the method builds no MSA for this target (Boltz, RNA only)
  !! AF3 msa_s and ColabFold-server msa_s are NOT comparable. Never plot them side by side. !!
  inference_s definition differs by method (column inference_s_definition):
    af3:    AF3's own "Running model inference with seed N took" (excludes container start, featurisation, I/O)
    boltz2, af2m: wall time of the per-seed run with the MSA cached (includes model load and I/O)
  run_wall_s: the whole per-seed Slurm step, for every method.

confidence.csv: long format: target, method, seed, sample, metric, value.
"""
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
T = ROOT / "results/tables"

cache = defaultdict(set)
for r in csv.DictReader(open(ROOT / "msas/msa_server_cache_status.tsv"), delimiter="\t"):
    cache[(r["target"], r["method"])].add(r["cache_status"])


def msa_source(t, m):
    if m == "af3":
        return "af3_local_jackhmmer"
    s = cache.get((t, m), set())
    if s == {"cache_hit"}:
        return "colabfold_server_cached"
    if s == {"fresh_search"}:
        return "colabfold_server_fresh"
    if s == {"no_msa(rna_only)"}:
        return "none"
    if {"cache_hit", "fresh_search"} <= s:
        return "colabfold_server_mixed"
    return "unknown"


INF_DEF = {"af3": "af3_model_inference_log", "boltz2": "per_seed_run_wall_msa_cached",
           "af2m": "per_seed_run_wall_msa_cached"}

rows = list(csv.DictReader(open(T / "manifest.csv")))
runs = {}
for r in rows:
    k = (r["target"], r["method"], r["seed"])
    runs.setdefault(k, r)
with open(T / "timing.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["target", "method", "seed", "msa_source", "msa_s", "msa_s_scope", "inference_s",
                "inference_s_definition", "run_wall_s", "run_rc", "gpu"])
    for (t, m, s), r in sorted(runs.items()):
        w.writerow([t, m, s, msa_source(t, m), r["msa_s"], "per_target_once", r["run_inference_s"],
                    INF_DEF[m], r["run_wall_s"], r["run_rc"], "L40S"])

CONF = ["ptm", "iptm", "ranking_score", "confidence_score", "plddt_mean", "complex_iplddt",
        "ligand_iptm", "protein_iptm", "has_clash", "fraction_disordered"]
n = 0
with open(T / "confidence.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["target", "method", "seed", "sample", "metric", "value"])
    for r in rows:
        for c in CONF:
            if r.get(c) not in (None, "", "None"):
                w.writerow([r["target"], r["method"], r["seed"], r["sample"], c, r[c]])
                n += 1
print(f"timing.csv: {len(runs)} runs; confidence.csv: {n} rows")
