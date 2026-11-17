"""Build results/tables/summary.csv from scores.csv + confidence.csv.

  envs/p3_tools/bin/python -I scripts/make_summary.py

One row per target x method x metric with:
  n, mean, median (over all 25 predictions)
  top_ranked   value of the method's own top-ranked prediction among its 25
               (AF3 / AF2-M: max ranking_score; Boltz-2: max confidence_score)
  best_of_25   best value among the 25 (max, or min for lower-is-better metrics)
  ranking_loss best_of_25 minus top_ranked, oriented so it is >= 0
               (higher-is-better: best - top; lower-is-better: top - best)
  seed_mean_min / seed_mean_max / seed_range / seed_means
               mean of the metric within each seed (5 samples), then min, max and range across
               the 5 seeds. The 5 samples within a seed share the trunk/MSA pass and are not
               independent, so the across-seed range is the honest spread.
  is_headline  the headline metric for the target type:
               protein-protein dockq; protein-NA dockq_protein_na_mean; RNA tm_score;
               protein-ligand success (BiSyRMSD < 2 A and PoseBusters-valid)
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
T = ROOT / "results/tables"
PP = ["dockq", "fnat", "irmsd", "usalign_tm_mm1"]
NA = ["dockq_protein_na_mean", "ics_protein_na_mean", "ilddt_protein_na_mean", "lddt", "tm_score"]
LIG = ["success", "bisyrmsd", "lddt_pli", "pb_valid", "plif_recovery"]
METRICS = {"T1206": PP, "H1215": PP,
           "M1211": NA + ["rna_only_tm_score", "rna_only_lddt"], "M1276": NA,
           "R1211": ["tm_score", "lddt"],
           "L1002": LIG, "L2001": LIG, "L4001": LIG + ["success_fraction_copies"]}
HEADLINE = {"T1206": "dockq", "H1215": "dockq", "M1211": "dockq_protein_na_mean",
            "M1276": "dockq_protein_na_mean", "R1211": "tm_score",
            "L1002": "success", "L2001": "success", "L4001": "success"}
LOWER = {"irmsd", "bisyrmsd"}
KEY = ["target", "method", "seed", "sample"]


def main():
    s = pd.read_csv(T / "scores.csv")
    c = pd.read_csv(T / "confidence.csv")
    rank = (c[c.metric.isin(["ranking_score", "confidence_score"])]
            .pivot_table(index=KEY, values="value", aggfunc="first")["value"].rename("rank_conf"))
    ttype = s.drop_duplicates("target").set_index("target")["target_type"]
    rows = []
    for t, mets in METRICS.items():
        if t not in set(s.target):
            continue
        w = s[(s.target == t) & s.metric.isin(mets)].pivot_table(index=KEY, columns="metric", values="value")
        w = w.join(rank)
        for m, g in w.groupby(level="method"):
            top = g.loc[g.rank_conf.idxmax()]
            for met in mets:
                if met not in g:
                    continue
                v = g[met]
                lower = met in LOWER
                best = v.min() if lower else v.max()
                seed_means = v.groupby(level="seed").mean()
                rows.append(dict(
                    target=t, target_type=ttype[t], method=m, metric=met,
                    is_headline=int(met == HEADLINE[t]), direction="lower_better" if lower else "higher_better",
                    n=int(v.notna().sum()), mean=v.mean(), median=v.median(),
                    top_ranked=top[met], best_of_25=best,
                    ranking_loss=(top[met] - best) if lower else (best - top[met]),
                    seed_mean_min=seed_means.min(), seed_mean_max=seed_means.max(),
                    seed_range=seed_means.max() - seed_means.min(),
                    seed_means=";".join(f"s{k}={x:.3f}" for k, x in seed_means.items())))
    d = pd.DataFrame(rows)
    num = d.select_dtypes("number").columns.drop(["n", "is_headline"])
    d[num] = d[num].round(4)
    d.to_csv(T / "summary.csv", index=False)
    print(f"summary.csv: {len(d)} rows")


if __name__ == "__main__":
    main()
