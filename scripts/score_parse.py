"""Stage D: parse raw scorer outputs into results/tables/scores.csv (long format).

  envs/p3_tools/bin/python -I scripts/score_parse.py

Columns: target, target_type, method, seed, sample, metric, value

Metrics
  every polymer target (OST 2.12.0 compare-structures; common residues only):
    lddt, ilddt, tm_score (OST built-in USalign), qs_global, ics, ips, rmsd (rigid),
    dockq_wave, dockq_ave, plus per-interface dockq / irmsd / lrmsd / fnat named
    "<metric>[<refchain1>-<refchain2>]" (reference chain ids)
  protein-protein (T1206, H1215):
    dockq (the single interface), irmsd, usalign_tm_mm1 (USalign -mm 1 -ter 0, normalised
    by the reference, i.e. TM2)
  protein-RNA / protein-DNA (M1211, M1276), HEADLINE = interface, not the whole complex:
    dockq_protein_na_mean: mean DockQ over protein-nucleic-acid interfaces only
    irmsd_protein_na_mean: likewise for interface RMSD
  RNA only (R1211): tm_score, lddt
  M1211 RNA chain scored alone (chain C vs native chain 0, same OST options as R1211):
    rna_only_tm_score, rna_only_lddt  -> directly comparable with R1211 tm_score / lddt
  protein-ligand (L1002, L2001; OST compare-ligand-structures + PoseBusters dock mode):
    bisyrmsd, lddt_pli, lddt_lp, bb_rmsd, pb_valid (1 = every PoseBusters check passed),
    pb_n_failed, success (1 = bisyrmsd < 2 A AND pb_valid)  <- HEADLINE
  failures: scoring_failed = 1 (raw error in predictions/scoring_work/out/<tag>.err);
            ligand_unassigned = 1 if OST could not assign the model ligand.
"""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
W = ROOT / "predictions/scoring_work"
OUT = W / "out"
TARGETS = json.load(open(ROOT / "data/targets/targets.json"))
NA_REF_CHAINS = {"M1211": {"0"}, "M1276": {"0"}}  # nucleic-acid chain ids in the natives
PP = {"T1206", "H1215"}
POLY_KEYS = ["lddt", "ilddt", "tm_score", "qs_global", "ics", "ips", "rmsd", "dockq_wave", "dockq_ave"]


def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def polymer_metrics(t, tag):
    d = json.load(open(OUT / f"{tag}.ost.json"))
    if d.get("status") != "SUCCESS":
        return {"scoring_failed": 1}
    m = {k: d.get(k) for k in POLY_KEYS}
    ints = d.get("dockq_interfaces") or []
    dq, ir, lr, fn = d.get("dockq") or [], d.get("irmsd") or [], d.get("lrmsd") or [], d.get("fnat") or []
    na_dq, na_ir = [], []
    for i, it in enumerate(ints):
        lab = f"[{it[0]}-{it[1]}]"  # reference chain pair
        m[f"dockq{lab}"], m[f"irmsd{lab}"] = dq[i], ir[i]
        m[f"lrmsd{lab}"], m[f"fnat{lab}"] = lr[i], fn[i]
        if t in NA_REF_CHAINS and (it[0] in NA_REF_CHAINS[t]) != (it[1] in NA_REF_CHAINS[t]):
            na_dq.append(dq[i])
            na_ir.append(ir[i])
    if t in PP:
        if len(ints) == 1:
            m["dockq"], m["irmsd"] = dq[0], ir[0]
        else:
            m["dockq"], m["irmsd"] = None, None  # should not happen for a dimer; left blank
        us = OUT / f"{tag}.usalign.txt"
        if us.exists():
            rows = [l.split("\t") for l in us.read_text().splitlines() if l and not l.startswith("#")]
            if rows:
                m["usalign_tm_mm1"] = float(rows[0][3])  # TM2: normalised by structure 2 = reference
    if t in NA_REF_CHAINS:
        m["dockq_protein_na_mean"] = mean(na_dq)
        m["irmsd_protein_na_mean"] = mean(na_ir)
        m["n_protein_na_interfaces"] = len(na_dq)
    return m


def ligand_metrics(tag):
    m = {}
    d = json.load(open(OUT / f"{tag}.ostlig.json"))
    if d.get("status") != "SUCCESS":
        return {"scoring_failed": 1}
    r = d.get("rmsd", {}).get("assigned_scores", [])
    lp = d.get("lddt_pli", {}).get("assigned_scores", [])
    if not r:
        m["ligand_unassigned"] = 1
    else:
        m["bisyrmsd"], m["lddt_lp"], m["bb_rmsd"] = r[0]["score"], r[0].get("lddt_lp"), r[0].get("bb_rmsd")
    if lp:
        m["lddt_pli"] = lp[0]["score"]
    pb = OUT / f"{tag}.posebusters.csv"
    rows = list(csv.DictReader(open(pb))) if pb.exists() else []
    if rows:
        checks = {k: v for k, v in rows[0].items() if v in ("True", "False")}
        failed = [k for k, v in checks.items() if v == "False"]
        m["pb_valid"] = int(not failed)
        m["pb_n_failed"] = len(failed)
        for k in failed:
            m[f"pb_fail:{k}"] = 1
    if "bisyrmsd" in m and "pb_valid" in m:
        m["success"] = int(m["bisyrmsd"] < 2.0 and m["pb_valid"] == 1)
    return m


def main():
    rows, n_fail = [], 0
    for t in csv.DictReader(open(W / "tasks.tsv"), delimiter="\t"):
        tag, tg = t["tag"], t["target"]
        if (OUT / f"{tag}.err").exists():
            met = {"scoring_failed": 1}
        elif t["kind"] == "polymer":
            met = polymer_metrics(tg, tag)
        elif t["kind"] == "ligand":
            met = ligand_metrics(tag)
        elif t["kind"] == "rna_only":
            d = json.load(open(OUT / f"{tag}.ost.json"))
            met = ({"rna_only_tm_score": d.get("tm_score"), "rna_only_lddt": d.get("lddt")}
                   if d.get("status") == "SUCCESS" else {"scoring_failed": 1})
        else:
            met = {"scoring_failed": 1}
        n_fail += met.get("scoring_failed", 0)
        for k, v in met.items():
            if v is not None:
                rows.append([tg, TARGETS[tg]["type"], t["method"], t["seed"], t["sample"], k,
                             round(v, 4) if isinstance(v, float) else v])
    out = ROOT / "results/tables/scores.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["target", "target_type", "method", "seed", "sample", "metric", "value"])
        w.writerows(rows)
    print(f"{len(rows)} rows -> {out}; {n_fail} predictions failed scoring")


if __name__ == "__main__":
    main()
