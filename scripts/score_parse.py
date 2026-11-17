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
  added in part 3 (2026-11-16):
    T1206/H1215: fnat (single interface, from the DockQ output)
    M1211/M1276: ics[<r1>-<r2>], ilddt[<r1>-<r2>] per protein-NA interface, and
      ics_protein_na_mean, ilddt_protein_na_mean (same restriction as dockq_protein_na_mean;
      ilddt from scripts/score_na_interfaces.py)
    ligands: plif_recovery, plif_n_native, plif_n_recovered (scripts/score_prolif.py)
    L4001 (2 ligand copies): per-copy bisyrmsd/lddt_pli/pb_valid/success[ligand0|1];
      bisyrmsd and lddt_pli = mean over assigned copies; pb_valid = all copies valid;
      success = 1 only if EVERY reference copy is placed < 2 A and PB-valid;
      success_fraction_copies; ligand_unassigned = number of unassigned reference copies
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


def polymer_metrics(t, tag):  # noqa: C901
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
            m["fnat"] = fn[0]  # 3a: fraction of native contacts, from the DockQ output
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
        # 3b: ICS and interface lDDT, restricted to protein-NA interfaces like dockq_protein_na_mean
        na_ics, na_il = [], []
        for (r1, r2), ics in zip(d.get("contact_reference_interfaces") or [], d.get("per_interface_ics") or []):
            if (r1 in NA_REF_CHAINS[t]) == (r2 in NA_REF_CHAINS[t]):
                continue
            m[f"ics[{r1}-{r2}]"] = ics
            na_ics.append(ics)
            f = OUT / f"{tag}.naiface_{r1}-{r2}.json"
            j = json.load(open(f)) if f.exists() else {}
            il = j.get("ilddt") if j.get("status") == "SUCCESS" else None
            m[f"ilddt[{r1}-{r2}]"] = il
            na_il.append(il)
        m["ics_protein_na_mean"] = mean(na_ics)
        m["ilddt_protein_na_mean"] = mean(na_il) if None not in na_il else None
    return m


def pb_rows(tag):
    pb = OUT / f"{tag}.posebusters.csv"
    rows = list(csv.DictReader(open(pb))) if pb.exists() else []
    out = {}
    for r in rows:
        checks = {k: v for k, v in r.items() if v in ("True", "False")}
        out[r.get("file")] = [k for k, v in checks.items() if v == "False"]
    return out


def prolif_metrics(tag):
    f = OUT / f"{tag}.prolif.json"
    if not f.exists():
        return {}
    d = json.load(open(f))
    return {"plif_recovery": d["recovery"], "plif_n_native": d["n_native"], "plif_n_recovered": d["n_recovered"]}


def ligand_metrics(tag, n_ref_ligands):
    m = {}
    d = json.load(open(OUT / f"{tag}.ostlig.json"))
    if d.get("status") != "SUCCESS":
        return {"scoring_failed": 1}
    r = d.get("rmsd", {}).get("assigned_scores", [])
    lp = d.get("lddt_pli", {}).get("assigned_scores", [])
    pbr = pb_rows(tag)
    if n_ref_ligands == 1:  # L1002, L2001: unchanged definitions
        if not r:
            m["ligand_unassigned"] = 1
        else:
            m["bisyrmsd"], m["lddt_lp"], m["bb_rmsd"] = r[0]["score"], r[0].get("lddt_lp"), r[0].get("bb_rmsd")
        if lp:
            m["lddt_pli"] = lp[0]["score"]
        if pbr:
            failed = next(iter(pbr.values()))
            m["pb_valid"] = int(not failed)
            m["pb_n_failed"] = len(failed)
            for k in failed:
                m[f"pb_fail:{k}"] = 1
        if "bisyrmsd" in m and "pb_valid" in m:
            m["success"] = int(m["bisyrmsd"] < 2.0 and m["pb_valid"] == 1)
    else:  # L4001: one row set per reference ligand copy + pooled values (see CLAUDE.md)
        lp_by_ref = {x["reference_ligand"]: x["score"] for x in lp}
        copy_ok, rmsds, plis, pbv = [], [], [], []
        for x in r:
            ref = Path(x["reference_ligand"]).stem.split("_")[-1]  # ligand0 / ligand1
            failed = pbr.get(x["model_ligand"], None)
            m[f"bisyrmsd[{ref}]"] = x["score"]
            m[f"lddt_pli[{ref}]"] = lp_by_ref.get(x["reference_ligand"])
            if failed is not None:
                m[f"pb_valid[{ref}]"] = int(not failed)
                for k in failed:
                    m[f"pb_fail:{k}[{ref}]"] = 1
            ok = int(x["score"] < 2.0 and failed is not None and not failed)
            m[f"success[{ref}]"] = ok
            copy_ok.append(ok)
            rmsds.append(x["score"])
            plis.append(lp_by_ref.get(x["reference_ligand"]))
            pbv.append(failed is not None and not failed)
        n_unassigned = n_ref_ligands - len(r)
        m["ligand_unassigned"] = n_unassigned
        copy_ok += [0] * n_unassigned
        m["bisyrmsd"] = mean(rmsds)          # mean over ASSIGNED copies
        m["lddt_pli"] = mean(plis)
        m["pb_valid"] = int(len(pbv) == n_ref_ligands and all(pbv))
        m["success"] = int(all(copy_ok))     # every reference copy placed < 2 A and PB-valid
        m["success_fraction_copies"] = sum(copy_ok) / n_ref_ligands
    m.update(prolif_metrics(tag))
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
            met = ligand_metrics(tag, len(t["reference_ligand"].split(";")))
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
