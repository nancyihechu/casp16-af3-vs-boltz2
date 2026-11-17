"""Figures for the README: every figure in a light and a dark version.

  envs/p3_tools/bin/python -I scripts/make_figures.py

Reads only results/tables/{scores,confidence,timing,msa_depth}.csv; writes results/figures/*.png.

Colour: methods use categorical slots 1-3 of the reference palette, validated with the dataviz
validator (all pairs, both modes: worst CVD dE 9.2 light / 9.4 dark, normal-vision 24.0 / 20.9).
Light aqua is < 3:1 on the light surface, so every method is also named in text and occupies its
own lane / marker shape (identity is never colour alone).
Quality thresholds used for the "fingerprint" figure (stated on the figure):
  DockQ >= 0.23 (CAPRI 'acceptable'), RNA TM-score >= 0.45, ligand success = BiSyRMSD < 2 A and
  PoseBusters-valid.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
T = ROOT / "results/tables"
FIG = ROOT / "results/figures"
FIG.mkdir(parents=True, exist_ok=True)

THEMES = {
    "light": dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", muted="#898781", grid="#e1e0d9",
                  axis="#c3c2b7", wash="#f0efec",
                  m={"af3": "#2a78d6", "boltz2": "#eb6834", "af2m": "#1baf7a"}),
    "dark": dict(surface="#1a1a19", ink="#f0efec", ink2="#c3c2b7", muted="#898781", grid="#2c2c2a",
                 axis="#383835", wash="#2a2a28",
                 m={"af3": "#3987e5", "boltz2": "#d95926", "af2m": "#199e70"}),
}
NAME = {"af3": "AlphaFold3", "boltz2": "Boltz-2", "af2m": "AF2-Multimer"}
MARK = {"af3": "o", "boltz2": "D", "af2m": "s"}
METHODS = ["af3", "boltz2", "af2m"]
TARGETS = ["T1206", "H1215", "M1211", "M1276", "R1211", "L1002", "L2001", "L4001"]
DESC = {"T1206": "protein homodimer", "H1215": "nanobody : mNeonGreen", "M1211": "protein : RNA",
        "M1276": "protein : ssDNA", "R1211": "RNA alone", "L1002": "chymase + ligand",
        "L2001": "cathepsin G + ligand", "L4001": "Mpro dimer + 2 ligands"}
HEAD = {"T1206": "dockq", "H1215": "dockq", "M1211": "dockq_protein_na_mean", "M1276": "dockq_protein_na_mean",
        "R1211": "tm_score", "L1002": "bisyrmsd", "L2001": "bisyrmsd", "L4001": "bisyrmsd"}
HEADLABEL = {"dockq": "DockQ", "dockq_protein_na_mean": "interface DockQ (protein–NA)", "tm_score": "TM-score",
             "bisyrmsd": "ligand BiSyRMSD (Å)"}
KEY = ["target", "method", "seed", "sample"]

plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Inter", "Helvetica Neue", "Arial", "DejaVu Sans"],
                     "font.size": 10, "axes.linewidth": 0.8, "svg.fonttype": "none"})


# ---------------------------------------------------------------- data
def load():
    s = pd.read_csv(T / "scores.csv")
    c = pd.read_csv(T / "confidence.csv")
    w = s.pivot_table(index=KEY, columns="metric", values="value")
    conf = c[c.metric.isin(["ranking_score", "confidence_score"])].pivot_table(index=KEY, values="value")["value"]
    w = w.join(conf.rename("conf")).reset_index()
    w["head"] = w.apply(lambda r: r[HEAD[r["target"]]], axis=1)

    def good(r):
        if r["target"].startswith("L"):
            return r["success"] == 1
        if r["target"] == "R1211":
            return r["tm_score"] >= 0.45
        return r["head"] >= 0.23

    w["good"] = w.apply(good, axis=1)
    # rank of each sample within its seed by the method's own confidence (0 = most confident)
    w["rank_in_seed"] = w.groupby(["target", "method", "seed"])["conf"].rank(ascending=False, method="first") - 1
    # exactly ONE top-ranked prediction per target x method, chosen exactly as summary.csv does
    # (first maximum in (seed, sample) order; AF3 per-sample ranking scores are stored to 2 decimals and tie)
    w = w.sort_values(KEY).reset_index(drop=True)
    top_idx = w.groupby(["target", "method"])["conf"].idxmax()
    w["is_top"] = False
    w.loc[top_idx.values, "is_top"] = True
    return w


def style(ax, th, grid_axis="x"):
    ax.set_facecolor(th["surface"])
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(th["axis"])
    ax.tick_params(colors=th["muted"], length=0, labelsize=8.5)
    ax.grid(axis=grid_axis, color=th["grid"], linewidth=0.7)
    ax.set_axisbelow(True)


def fig_new(th, w, h):
    f = plt.figure(figsize=(w, h), facecolor=th["surface"])
    return f


def title(f, th, t, sub, y=0.975):
    f.text(0.035, y, t, fontsize=15, weight="bold", color=th["ink"], va="top")
    f.text(0.035, y - 0.045 * (8.5 / f.get_figheight()), sub, fontsize=9.5, color=th["ink2"], va="top")


def save(f, name, mode):
    p = FIG / f"{name}_{mode}.png"
    f.savefig(p, dpi=200, facecolor=f.get_facecolor())
    plt.close(f)
    return p


# ---------------------------------------------------------------- 1. fingerprint (hero)
def fig_fingerprint(w, th, mode):
    f = fig_new(th, 11, 9.2)
    title(f, th, "Every prediction, one tile",
          "Each block is one method on one target: 25 predictions = 5 seeds (rows) × 5 samples, ordered left→right by "
          "the method's own confidence.\nFilled = meets the quality bar   ·   DockQ ≥ 0.23 (CAPRI acceptable) · RNA "
          "TM-score ≥ 0.45 · ligand BiSyRMSD < 2 Å and PoseBusters-valid", y=0.975)
    ax = f.add_axes([0.03, 0.02, 0.94, 0.86])
    ax.set_facecolor(th["surface"])
    ax.axis("off")
    cell, gap, block_gap = 1.0, 0.18, 1.6
    bw = 5 * cell + 4 * gap
    col_x = {m: 12.5 + i * (bw + block_gap + 3.2) for i, m in enumerate(METHODS)}
    row_h = bw + 1.15
    groups = [("protein–protein", ["T1206", "H1215"]), ("protein–nucleic acid", ["M1211", "M1276"]),
              ("RNA", ["R1211"]), ("protein–ligand", ["L1002", "L2001", "L4001"])]
    y = 0.0
    for m in METHODS:
        ax.text(col_x[m] + bw / 2, 1.1, NAME[m], ha="center", va="bottom", fontsize=11.5, weight="bold",
                color=th["ink"])
        ax.plot([col_x[m], col_x[m] + bw], [0.75, 0.75], color=th["m"][m], lw=3, solid_capstyle="round")
    y = -0.6
    for gname, tg in groups:
        ax.text(0, y + 0.1, gname.upper(), fontsize=8, color=th["muted"], va="top", weight="bold")
        y -= 0.95
        for t in tg:
            ytop = y
            ax.text(0, ytop - bw / 2 + 0.45, t, fontsize=12, weight="bold", color=th["ink"], va="center")
            ax.text(0, ytop - bw / 2 - 0.55, DESC[t], fontsize=8.5, color=th["ink2"], va="center")
            for m in METHODS:
                g = w[(w.target == t) & (w.method == m)]
                x0 = col_x[m]
                if g.empty:
                    bg = FancyBboxPatch((x0, ytop - bw), bw, bw, boxstyle="round,pad=0,rounding_size=0.35",
                                        fc=th["wash"], ec="none")
                    ax.add_patch(bg)
                    ax.text(x0 + bw / 2, ytop - bw / 2, "not run\n(protein-only\nmethod)", ha="center", va="center",
                            fontsize=7.5, color=th["muted"])
                    continue
                for r in g.itertuples():
                    cx = x0 + r.rank_in_seed * (cell + gap)
                    cy = ytop - (r.seed - 1) * (cell + gap) - cell
                    if r.good:
                        p = FancyBboxPatch((cx, cy), cell, cell, boxstyle="round,pad=0,rounding_size=0.18",
                                           fc=th["m"][m], ec="none")
                    else:
                        p = FancyBboxPatch((cx + 0.06, cy + 0.06), cell - 0.12, cell - 0.12,
                                           boxstyle="round,pad=0,rounding_size=0.16", fc="none",
                                           ec=th["axis"], lw=1.0)
                    ax.add_patch(p)
                    if r.is_top:  # the method's single top-ranked prediction
                        ax.plot(cx + cell / 2, cy + cell / 2, marker="o", ms=3.2,
                                color=th["surface"] if r.good else th["ink2"], zorder=5)
                n = int(g.good.sum())
                ax.text(x0 + bw + 0.45, ytop - bw / 2, f"{n}", fontsize=13, weight="bold", color=th["ink"],
                        va="center", ha="left")
                ax.text(x0 + bw + 0.45, ytop - bw / 2 - 0.85, "/25", fontsize=8, color=th["muted"], va="center")
            y -= row_h
        y -= 0.35
    # legend
    ly = y + 0.2
    lx = 0
    ax.add_patch(FancyBboxPatch((lx, ly - 0.8), 0.8, 0.8, boxstyle="round,pad=0,rounding_size=0.15",
                                fc=th["ink2"], ec="none"))
    ax.text(lx + 1.2, ly - 0.4, "meets the bar", fontsize=8.5, color=th["ink2"], va="center")
    ax.add_patch(FancyBboxPatch((lx + 9, ly - 0.75), 0.7, 0.7, boxstyle="round,pad=0,rounding_size=0.15",
                                fc="none", ec=th["axis"], lw=1.0))
    ax.text(lx + 10.1, ly - 0.4, "does not", fontsize=8.5, color=th["ink2"], va="center")
    ax.plot(lx + 16.4, ly - 0.4, marker="o", ms=4, color=th["ink2"])
    ax.text(lx + 17.0, ly - 0.4, "the method's own top-ranked prediction (one per block)", fontsize=8.5,
            color=th["ink2"], va="center")
    ax.set_xlim(-0.2, col_x["af2m"] + bw + 2.4)
    ax.set_ylim(ly - 1.4, 2.4)
    ax.set_aspect("equal")
    return save(f, "01_fingerprint", mode)


# ---------------------------------------------------------------- 2. distributions
def fig_strips(w, th, mode):
    f = fig_new(th, 12, 10)
    title(f, th, "All 25 predictions per method, against the native",
          "Each dot is one prediction (seeds stacked within a lane). Ring = the method's own top-ranked pick; solid bar = "
          "mean of 25; grey rule = quality bar.")
    rng = np.random.default_rng(7)
    for i, t in enumerate(TARGETS):
        ax = f.add_axes([0.10 + (i % 2) * 0.48, 0.73 - (i // 2) * 0.215, 0.37, 0.13])
        style(ax, th)
        met = HEAD[t]
        meths = [m for m in METHODS if not w[(w.target == t) & (w.method == m)].empty]
        for j, m in enumerate(meths):
            g = w[(w.target == t) & (w.method == m)]
            yy = len(meths) - 1 - j
            jit = (g.seed.values - 3) * 0.075 + rng.uniform(-0.03, 0.03, len(g))
            ax.scatter(g["head"], yy + jit, s=20, color=th["m"][m], marker=MARK[m], lw=0, alpha=0.9, zorder=3)
            top = g[g.is_top].iloc[0]
            ax.scatter([top["head"]], [yy], s=110, facecolor="none", edgecolor=th["ink"], lw=1.4, zorder=4)
            mu = g["head"].mean()
            ax.plot([mu, mu], [yy - 0.33, yy + 0.33], color=th["ink"], lw=2, solid_capstyle="round", zorder=4)
            ax.text(-0.02, yy, NAME[m], transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=8.5, color=th["ink2"])
        ax.set_yticks([])
        ax.set_ylim(-0.6, len(meths) - 0.4)
        if met == "bisyrmsd":
            ax.set_xlim(0.3, 5.0)
            ax.axvline(2.0, color=th["ink2"], lw=1, zorder=2)
            ax.text(2.0, len(meths) - 0.35, " 2 Å", color=th["ink2"], fontsize=8, va="bottom")
            ax.text(1.0, 1.18, "← better", transform=ax.transAxes, ha="right", fontsize=8, color=th["muted"])
        else:
            ax.set_xlim(0, 1)
            bar = 0.45 if t == "R1211" else 0.23
            ax.axvline(bar, color=th["ink2"], lw=1, zorder=2)
            ax.text(bar, len(meths) - 0.35, f" {bar}", color=th["ink2"], fontsize=8, va="bottom")
            ax.text(1.0, 1.18, "better →", transform=ax.transAxes, ha="right", fontsize=8, color=th["muted"])
        ax.text(0, 1.18, t, transform=ax.transAxes, fontsize=12, weight="bold", color=th["ink"])
        ax.text(0.17, 1.18, f"{DESC[t]}  ·  {HEADLABEL[met]}", transform=ax.transAxes, fontsize=8.8,
                color=th["ink2"], va="baseline")
    return save(f, "02_distributions", mode)


# ---------------------------------------------------------------- 3. ranking anatomy & seed spread
def fig_anatomy(w, th, mode):
    f = fig_new(th, 12, 7.8)
    title(f, th, "How much the ranking and the seed matter",
          "Per method and target: the band spans the five per-seed means (honest spread), the bar is the mean of 25, "
          "◆ is the method's top-ranked pick, ▲ the best of 25.\nThe gap ◆→▲ is the ranking loss.")
    panels = [("DockQ  (protein–protein, protein–NA interface)", ["T1206", "H1215", "M1211", "M1276"], (0, 1)),
              ("TM-score  (RNA)", ["R1211"], (0, 1)),
              ("BiSyRMSD, Å  (ligands; lower is better)", ["L1002", "L2001", "L4001"], (0.3, 3.0))]
    heights = [len(p[1]) for p in panels]
    tops = [0.80, 0.41, 0.25]
    hs = [0.32, 0.075, 0.20]
    for (lab, tg, xl), top, h in zip(panels, tops, hs):
        ax = f.add_axes([0.20, top - h, 0.75, h])
        style(ax, th)
        rows = [(t, m) for t in tg for m in METHODS if not w[(w.target == t) & (w.method == m)].empty]
        for k, (t, m) in enumerate(rows):
            g = w[(w.target == t) & (w.method == m)]
            y = len(rows) - 1 - k
            sm = g.groupby("seed")["head"].mean()
            ax.plot([sm.min(), sm.max()], [y, y], color=th["m"][m], alpha=0.30, lw=9, solid_capstyle="butt",
                    zorder=2)
            ax.scatter(sm.values, [y] * 5, s=9, color=th["m"][m], zorder=3, lw=0)
            mu = g["head"].mean()
            ax.plot([mu, mu], [y - 0.3, y + 0.3], color=th["ink"], lw=2, zorder=4, solid_capstyle="round")
            top_v = g[g.is_top]["head"].iloc[0]
            best = g["head"].min() if xl[0] > 0 else g["head"].max()
            ax.plot([top_v, best], [y, y], color=th["ink2"], lw=1, zorder=3)
            ax.scatter([top_v], [y], marker="D", s=36, color=th["ink"], zorder=5, edgecolor=th["surface"], lw=1.5)
            ax.scatter([best], [y], marker="^", s=46, color=th["ink"], zorder=5, edgecolor=th["surface"], lw=1.5)
            ax.text(-0.03, y, f"{t}  {NAME[m]}", transform=ax.get_yaxis_transform(), ha="right", va="center",
                    fontsize=8.5, color=th["ink"] if k == 0 or rows[k - 1][0] != t else th["ink2"])
            ax.plot([-0.022, -0.01], [y, y], transform=ax.get_yaxis_transform(), color=th["m"][m], lw=4,
                    clip_on=False, solid_capstyle="round")
        ax.set_yticks([])
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_xlim(*xl)
        if xl[0] > 0:
            ax.axvline(2.0, color=th["ink2"], lw=1, zorder=1)
        ax.text(0, 1.02, lab, transform=ax.transAxes, fontsize=9.5, weight="bold", color=th["ink2"], va="bottom")
    return save(f, "03_ranking_and_seed_spread", mode)


# ---------------------------------------------------------------- 4. confidence vs accuracy
def fig_confidence(w, th, mode):
    f = fig_new(th, 12, 7.6)
    title(f, th, "Did the models know when they were right?",
          "x = the method's own confidence (AF3 / AF2-Multimer ranking score, Boltz-2 confidence score), "
          "y = accuracy against the native. One point per prediction.")
    for i, t in enumerate(TARGETS):
        ax = f.add_axes([0.075 + (i % 4) * 0.235, 0.49 - (i // 4) * 0.39, 0.175, 0.27])
        style(ax, th, grid_axis="both")
        for m in METHODS:
            g = w[(w.target == t) & (w.method == m)]
            if g.empty:
                continue
            ax.scatter(g.conf, g["head"], s=22, marker=MARK[m], color=th["m"][m], alpha=0.85, lw=0.6,
                       edgecolor=th["surface"], zorder=3)
        ax.set_xlim(0, 1)
        met = HEAD[t]
        if met == "bisyrmsd":
            ax.set_ylim(5.0, 0.3)  # inverted: up = better
            ax.axhline(2.0, color=th["ink2"], lw=0.9)
        else:
            ax.set_ylim(0, 1)
            ax.axhline(0.45 if t == "R1211" else 0.23, color=th["ink2"], lw=0.9)
        ax.text(0, 1.06, t, transform=ax.transAxes, fontsize=11, weight="bold", color=th["ink"])
        short = {"dockq": "DockQ", "dockq_protein_na_mean": "interface DockQ", "tm_score": "TM-score",
                 "bisyrmsd": "BiSyRMSD, Å (inverted)"}[met]
        ax.set_ylabel(short, color=th["muted"], fontsize=8)
    f.text(0.5, 0.025, "method's own confidence  →", ha="center", fontsize=9, color=th["muted"])
    handles = [Line2D([], [], marker=MARK[m], ls="", color=th["m"][m], ms=7, label=NAME[m]) for m in METHODS]
    leg = f.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.97, 0.90), ncol=3, frameon=False,
                   fontsize=9.5, labelcolor=th["ink2"])
    return save(f, "04_confidence_vs_accuracy", mode)


# ---------------------------------------------------------------- 5. RNA alone vs in complex
def fig_rna(w, th, mode):
    f = fig_new(th, 10, 5.2)
    title(f, th, "The same 90-nt RNA, alone and inside the protein complex",
          "R1211 (RNA alone) vs the RNA chain of M1211 scored on its own against identical native coordinates. "
          "Dots = 25 predictions; line joins the means.")
    for i, (met_a, met_c, lab) in enumerate([("tm_score", "rna_only_tm_score", "TM-score"),
                                             ("lddt", "rna_only_lddt", "lDDT")]):
        ax = f.add_axes([0.08 + i * 0.47, 0.12, 0.38, 0.62])
        style(ax, th, grid_axis="y")
        rng = np.random.default_rng(3)
        for k, m in enumerate(["af3", "boltz2"]):
            a = w[(w.target == "R1211") & (w.method == m)][met_a]
            c = w[(w.target == "M1211") & (w.method == m)][met_c]
            off = -0.09 if k == 0 else 0.09
            ax.scatter(0 + off + rng.uniform(-0.04, 0.04, len(a)), a, s=16, color=th["m"][m], marker=MARK[m],
                       alpha=0.55, lw=0)
            ax.scatter(1 + off + rng.uniform(-0.04, 0.04, len(c)), c, s=16, color=th["m"][m], marker=MARK[m],
                       alpha=0.55, lw=0)
            ax.plot([0 + off, 1 + off], [a.mean(), c.mean()], color=th["m"][m], lw=2.2, zorder=4)
            for x, v in [(0 + off, a.mean()), (1 + off, c.mean())]:
                ax.scatter([x], [v], s=60, color=th["m"][m], marker=MARK[m], edgecolor=th["surface"], lw=2,
                           zorder=5)
            ax.plot([0.02], [0.955 - 0.07 * k], transform=ax.transAxes, marker=MARK[m], color=th["m"][m], ms=7,
                    ls="")
            ax.text(0.06, 0.955 - 0.07 * k, f"{NAME[m]}   mean {a.mean():.3f} → {c.mean():.3f}",
                    transform=ax.transAxes, va="center", fontsize=8.8, color=th["ink"])
        ax.set_xticks([0, 1], ["alone\n(R1211)", "in complex\n(M1211)"], color=th["ink2"], fontsize=9)
        ax.set_xlim(-0.4, 1.4)
        ax.set_ylim(0.3, 0.92)
        ax.set_yticks(np.round(np.arange(0.3, 0.81, 0.1), 1))
        ax.text(0, 1.04, lab, transform=ax.transAxes, fontsize=11, weight="bold", color=th["ink"])
    return save(f, "05_rna_alone_vs_complex", mode)


# ---------------------------------------------------------------- 6. ligands
def fig_ligands(w, th, mode):
    f = fig_new(th, 12, 4.8)
    title(f, th, "Ligand poses: geometry and interactions",
          "x = BiSyRMSD to the crystal pose (CASP16 metric), y = fraction of native protein–ligand interactions "
          "recovered (ProLIF).\nShaded = BiSyRMSD < 2 Å (success also requires PoseBusters validity; for L4001, both ligand copies).")
    for i, t in enumerate(["L1002", "L2001", "L4001"]):
        ax = f.add_axes([0.06 + i * 0.315, 0.13, 0.26, 0.58])
        style(ax, th, grid_axis="both")
        ax.axvspan(0, 2.0, color=th["wash"], zorder=0)
        for m in ["af3", "boltz2"]:
            g = w[(w.target == t) & (w.method == m)]
            bad = g[g.pb_valid < 1]
            ax.scatter(g.bisyrmsd, g.plif_recovery, s=26, marker=MARK[m], color=th["m"][m], alpha=0.85,
                       edgecolor=th["surface"], lw=0.6, zorder=3)
            if len(bad):
                ax.scatter(bad.bisyrmsd, bad.plif_recovery, s=120, facecolor="none", edgecolor=th["ink"], lw=1.2,
                           zorder=4)
                for r in bad.itertuples():
                    ax.annotate("PoseBusters\nfailure", (r.bisyrmsd, r.plif_recovery), xytext=(12, -22),
                                textcoords="offset points", fontsize=7.5, color=th["ink2"],
                                arrowprops=dict(arrowstyle="-", color=th["ink2"], lw=0.8))
        ax.set_xlim(0.3, 5.0)
        ax.set_ylim(0.0, 1.0)
        ax.text(0, 1.05, f"{t}", transform=ax.transAxes, fontsize=11.5, weight="bold", color=th["ink"])
        ax.text(0.2, 1.05, DESC[t], transform=ax.transAxes, fontsize=8.8, color=th["ink2"])
        ax.set_xlabel("BiSyRMSD (Å)", color=th["muted"], fontsize=8.5)
        if i == 0:
            ax.set_ylabel("interaction recovery", color=th["muted"], fontsize=8.5)
    handles = [Line2D([], [], marker=MARK[m], ls="", color=th["m"][m], ms=7, label=NAME[m]) for m in ["af3", "boltz2"]]
    f.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.97, 0.84), ncol=2, frameon=False,
             fontsize=9.5, labelcolor=th["ink2"])
    return save(f, "06_ligand_poses", mode)


# ---------------------------------------------------------------- 7. time
def fig_time(th, mode):
    t = pd.read_csv(T / "timing.csv")
    f = fig_new(th, 12, 5.4)
    title(f, th, "Where the time goes",
          "Left: AF3's MSA search, run once per target on CPU (local jackhmmer/nhmmer).\n"
          "Right: wall-clock of one seed on an NVIDIA L40S (5 samples, MSA cached, steady state; AF2-Multimer = job time ÷ 4 seeds).\n"
          "Boltz-2 and AF2-Multimer MSA times are deliberately not shown: a public server returned cached results, so they are not comparable.")
    msa = t[t.method == "af3"].drop_duplicates("target").set_index("target")["msa_s"] / 60
    st = t[t.stage.isin(["C_steady_state", "D_added_2026-11-16"])]
    st = st.assign(per_seed=st.run_wall_s / st.seeds_in_job)
    inf = st.groupby(["target", "method"])["per_seed"].mean()
    order = TARGETS
    ax = f.add_axes([0.10, 0.10, 0.32, 0.58])
    style(ax, th)
    ys = np.arange(len(order))[::-1]
    ax.barh(ys, [msa[x] for x in order], height=0.55, color=th["m"]["af3"], zorder=3)
    for y, x in zip(ys, order):
        ax.text(msa[x] + 0.4, y, f"{msa[x]:.0f} min", va="center", fontsize=8.5, color=th["ink"])
    ax.set_yticks(ys, order, color=th["ink"], fontsize=9.5)
    ax.set_xlim(0, 32)
    ax.text(0, 1.03, "AF3 MSA search, minutes (once per target)", transform=ax.transAxes, fontsize=9.5,
            weight="bold", color=th["ink2"])
    ax2 = f.add_axes([0.55, 0.10, 0.40, 0.58])
    style(ax2, th)
    for m in METHODS:
        xs, yv = [], []
        for y, x in zip(ys, order):
            if (x, m) in inf:
                xs.append(inf[(x, m)])
                yv.append(y + {"af3": 0.18, "boltz2": 0, "af2m": -0.18}[m])
        ax2.scatter(xs, yv, s=48, marker=MARK[m], color=th["m"][m], edgecolor=th["surface"], lw=1.2, zorder=3,
                    label=NAME[m])
    ax2.set_yticks(ys, order, color=th["ink"], fontsize=9.5)
    ax2.set_xlim(0, 340)
    ax2.text(0, 1.03, "Wall-clock per seed on GPU, seconds", transform=ax2.transAxes, fontsize=9.5,
             weight="bold", color=th["ink2"])
    ax2.legend(loc="lower right", frameon=False, fontsize=9, labelcolor=th["ink2"])
    return save(f, "07_time", mode)


def main():
    w = load()
    for mode, th in THEMES.items():
        for fn in (fig_fingerprint, fig_strips, fig_anatomy, fig_confidence, fig_rna, fig_ligands):
            print(fn(w, th, mode).name)
        print(fig_time(th, mode).name)


if __name__ == "__main__":
    main()
