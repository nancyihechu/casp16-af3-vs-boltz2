"""Emit the numeric tables for results/PART4_REPORT.md from the CSVs (no hand transcription).

  envs/p3_tools/bin/python -I scripts/make_report_tables.py > predictions/report_tables.md
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
T = ROOT / "results/tables"
s = pd.read_csv(T / "summary.csv")
M = {"af3": "AF3", "boltz2": "Boltz-2", "af2m": "AF2-Multimer"}


def f(x, nd=3):
    return "–" if pd.isna(x) else f"{x:.{nd}f}"


def table(targets, metrics):
    lines = ["| Target | Metric | Method | Mean of 25 | Top-ranked | Best of 25 | Ranking loss | Per-seed mean range (min–max) |",
             "|---|---|---|---|---|---|---|---|"]
    for t in targets:
        for met in metrics:
            for m in ["af3", "boltz2", "af2m"]:
                r = s[(s.target == t) & (s.metric == met) & (s.method == m)]
                if r.empty:
                    continue
                r = r.iloc[0]
                head = " (headline)" if r.is_headline else ""
                lines.append(f"| {t} | {met}{head} | {M[m]} | {f(r['mean'])} | {f(r.top_ranked)} | {f(r.best_of_25)} | "
                             f"{f(r.ranking_loss)} | {f(r.seed_mean_min)}–{f(r.seed_mean_max)} (range {f(r.seed_range)}) |")
    return "\n".join(lines)


print("### Protein–protein\n")
print(table(["T1206", "H1215"], ["dockq", "fnat", "irmsd", "usalign_tm_mm1"]))
print("\n### Protein–nucleic acid (interface metrics restricted to protein–NA interfaces)\n")
print(table(["M1211", "M1276"], ["dockq_protein_na_mean", "ics_protein_na_mean", "ilddt_protein_na_mean", "lddt", "tm_score"]))
print("\n### RNA (R1211 alone; M1211 RNA chain scored alone)\n")
print(table(["R1211"], ["tm_score", "lddt"]))
print(table(["M1211"], ["rna_only_tm_score", "rna_only_lddt"]).split("\n", 2)[2])
print("\n### Protein–ligand\n")
print(table(["L1002", "L2001", "L4001"], ["success", "bisyrmsd", "lddt_pli", "pb_valid", "plif_recovery", "success_fraction_copies"]))

d = pd.read_csv(T / "diagnostics/boltz2_potentials.csv")
print("\n### Diagnostic 1b: Boltz-2 potentials (L1002, L2001; seed 1; 5 samples each)\n")
print("| Target | Run | BiSyRMSD per sample (0–4) | Mean BiSyRMSD | PoseBusters-valid | Failed PB checks | Success | Run wall s |")
print("|---|---|---|---|---|---|---|---|")
for (t, lab), g in d.groupby(["target", "label"], sort=True):
    g = g.sort_values("sample")
    fails = sorted({c for x in g.pb_failed_checks.dropna() for c in str(x).split(";") if c})
    print(f"| {t} | {lab} | {', '.join(f(x) for x in g.bisyrmsd)} | {f(g.bisyrmsd.mean())} | "
          f"{int(g.pb_valid.sum())}/{len(g)} | {', '.join(fails) or 'none'} | {int(g.success.sum())}/{len(g)} | "
          f"{f(g.wall_s.iloc[0], 1)} |")
