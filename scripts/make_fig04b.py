"""Figure 04b (ipTM vs DockQ) and results/tables/confidence_vs_accuracy.csv. New files only.

  envs/p3_tools/bin/python -I scripts/make_fig04b.py

x = ipTM, y = DockQ (T1206, H1215) or dockq_protein_na_mean (M1211, M1276), values from scores.csv.
ipTM definition:
  T1206, H1215, M1276: the global ipTM each method reports (`iptm` in manifest.csv).
  M1211: protein-RNA chain-pair ipTM only, averaged over the protein-RNA pairs:
    AF3     chain_pair_iptm[A][C] and [B][C] (AF3's matrix is symmetric), mean of 2;
    Boltz-2 pair_chains_iptm 0->2, 2->0, 1->2, 2->1 (asymmetric: both directions), mean of 4.
Precision: AF3 and AF2-Multimer report ipTM to 2 decimals (ties possible); Boltz-2 full precision.
Top-ranked = the method's own ranking (AF3 full-precision ranking_score, AF2-M ranking_score,
Boltz-2 confidence_score), as everywhere else in this repo.
"Highest-ipTM prediction": if several predictions tie at the maximum ipTM, the reported DockQ is the
mean over the tied set and n_tied_at_max_iptm says how many.
Spearman: Pearson correlation of average ranks (ties get average ranks).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import make_figures as F  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

TARGETS = ["T1206", "H1215", "M1211", "M1276"]
YLAB = {"T1206": "DockQ", "H1215": "DockQ", "M1211": "interface DockQ (protein–RNA)",
        "M1276": "interface DockQ (protein–DNA)"}
THEMES = {
    "dark": F.THEMES["dark"],
    "light": dict(F.THEMES["light"], surface="#fcfcfb", ink="#3e4a59", ink2="#3e4a59", muted="#3e4a59"),
}


def m1211_pna_iptm(method, cpi):
    if method == "af3":
        mat = json.loads(cpi)  # rows/cols in chain order A, B, C(RNA)
        return np.mean([mat[0][2], mat[1][2]])
    d = json.loads(cpi)  # Boltz-2: {"0": {"0":..,"1":..,"2":..}, ...}; 0,1 protein, 2 RNA
    return np.mean([d["0"]["2"], d["2"]["0"], d["1"]["2"], d["2"]["1"]])


def data():
    w = F.load()
    man = pd.read_csv(ROOT / "results/tables/manifest.csv")[F.KEY + ["iptm", "chain_pair_iptm"]]
    w = w.drop(columns=[c for c in ["iptm"] if c in w.columns]).merge(man, on=F.KEY, how="left")
    w = w[w.target.isin(TARGETS)].copy()
    w["x_iptm"] = [m1211_pna_iptm(r.method, r.chain_pair_iptm) if r.target == "M1211" else r.iptm
                   for r in w.itertuples()]
    return w


def spearman(a, b):
    a, b = pd.Series(a).rank(), pd.Series(b).rank()
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def table(w):
    rows = []
    for t in TARGETS:
        for m in F.METHODS:
            g = w[(w.target == t) & (w.method == m)]
            if g.empty:
                continue
            mx = g.x_iptm.max()
            tied = g[np.isclose(g.x_iptm, mx)]
            rows.append(dict(
                target=t, method=m, n=len(g),
                iptm_definition=("protein-RNA chain-pair ipTM, mean over pairs" if t == "M1211" else "global ipTM"),
                iptm_precision=("full" if m == "boltz2" else "2 decimals"),
                mean_iptm=round(g.x_iptm.mean(), 4),
                spearman_iptm_vs_dockq=round(spearman(g.x_iptm, g["head"]), 4),
                spearman_ranking_score_vs_dockq=round(spearman(g.conf, g["head"]), 4),
                ranking_score_used=("confidence_score" if m == "boltz2" else "ranking_score"),
                dockq_top_ranked=round(float(g[g.is_top]["head"].iloc[0]), 4),
                dockq_at_max_iptm=round(float(tied["head"].mean()), 4),
                n_tied_at_max_iptm=len(tied)))
    d = pd.DataFrame(rows)
    d.to_csv(ROOT / "results/tables/confidence_vs_accuracy.csv", index=False)
    return d


def figure(w, th, mode):
    f = F.fig_new(th, 12, 7.6)
    F.title(f, th, "ipTM against interface accuracy",
            "x = ipTM; y = DockQ (T1206, H1215) or interface DockQ over protein–nucleic-acid interfaces (M1211, M1276).\n"
            "25 predictions per method; ring = the method's own top-ranked pick. ipTM = global ipTM, except M1211 = protein–RNA\n"
            "chain-pair ipTM averaged over pairs (AF3: A–C, B–C; Boltz-2: both directions, 4 values). AF3 / AF2-Multimer report ipTM to 2 decimals.")
    for i, t in enumerate(TARGETS):
        ax = f.add_axes([0.075 + (i % 2) * 0.47, 0.47 - (i // 2) * 0.39, 0.39, 0.27])
        F.style(ax, th, grid_axis="both")
        for m in F.METHODS:
            g = w[(w.target == t) & (w.method == m)]
            if g.empty:
                continue
            ax.scatter(g.x_iptm, g["head"], s=26, marker=F.MARK[m], color=th["m"][m], alpha=0.85, lw=0.6,
                       edgecolor=th["surface"], zorder=3)
            top = g[g.is_top].iloc[0]
            ax.scatter([top.x_iptm], [top["head"]], s=150, facecolor="none", edgecolor=th["ink"], lw=1.4, zorder=4)
        ax.axhline(0.23, color=th["ink2"], lw=0.9, zorder=2)
        ax.text(1.0, 0.235, "DockQ 0.23 ", ha="right", va="bottom", fontsize=7.5, color=th["ink2"],
                transform=ax.get_yaxis_transform())
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_ylabel(YLAB[t], color=th["muted"], fontsize=8.5)
        ax.text(0, 1.06, t, transform=ax.transAxes, fontsize=11.5, weight="bold", color=th["ink"])
        ax.text(0.17, 1.06, F.DESC[t], transform=ax.transAxes, fontsize=9, color=th["ink2"])
        ax.tick_params(colors=th["muted"])
    f.text(0.5, 0.025, "ipTM  →", ha="center", fontsize=9.5, color=th["muted"])
    handles = [Line2D([], [], marker=F.MARK[m], ls="", color=th["m"][m], ms=7, label=F.NAME[m]) for m in F.METHODS]
    handles.append(Line2D([], [], marker="o", ls="", markerfacecolor="none", markeredgecolor=th["ink"], ms=10,
                          label="top-ranked"))
    f.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.97, 0.84), ncol=4, frameon=False,
             fontsize=9.5, labelcolor=th["ink2"])
    return F.save(f, "04b_iptm_vs_dockq", mode)


def main():
    w = data()
    d = table(w)
    for mode, th in THEMES.items():
        print(figure(w, th, mode).name)
    pd.set_option("display.width", 250)
    print(d.drop(columns=["iptm_definition", "ranking_score_used"]).to_string(index=False))


if __name__ == "__main__":
    main()
