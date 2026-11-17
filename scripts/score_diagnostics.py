"""Part 1b: score the DIAGNOSTIC Boltz-2 potentials runs exactly like the benchmark ligand poses.

  envs/p3_tools/bin/python -I scripts/score_diagnostics.py

Reads predictions/diagnostics/boltz2_<label>/<T>/s<seed>/ (from scripts/run_boltz2_diag.slurm).
Uses the SAME functions as the benchmark pipeline (score_prep.coords_to_sdf / write_pdb,
score_tasks OST + PoseBusters commands, score_parse.ligand_metrics) but writes to
predictions/scoring_work/diag/ and results/tables/diagnostics/boltz2_potentials.csv only.
Nothing here touches manifest.csv, scores.csv or summary.csv.

Rows are labelled:
  benchmark_potentials_off  the benchmark Boltz-2 seed-1 poses (from scores.csv / confidence.csv)
  potentials_off_rerun      no --use_potentials (identical configuration; determinism check)
  potentials_on             --use_potentials
"""
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import score_parse  # noqa: E402
import score_prep  # noqa: E402
import score_tasks  # noqa: E402

DIAG = ROOT / "predictions/diagnostics"
WD = ROOT / "predictions/scoring_work/diag"
OUT = WD / "out"
OUT.mkdir(parents=True, exist_ok=True)
TARGETS = json.load(open(ROOT / "data/targets/targets.json"))
REF = {t: (ROOT / f"predictions/scoring_work/ref_{t}_ligand.sdf",
           next((ROOT / "predictions/scoring_work").glob(f"ref_{t}_receptor_chain*.pdb")))
       for t in ["L1002", "L2001"]}


def main():
    rows = []
    for d in sorted(DIAG.glob("boltz2_*/*/s*")):
        label, t, seed = d.parent.parent.name.replace("boltz2_", ""), d.parent.name, int(d.name[1:])
        timing = list(csv.DictReader(open(d / "timing.tsv"), delimiter="\t"))[0]
        pdir = d / f"boltz_results_{t}" / "predictions" / t
        smi = next(c["smiles"] for c in TARGETS[t]["chains"] if c["type"] == "ligand")
        rsdf, rpdb = REF[t]
        for cif in sorted(pdir.glob(f"{t}_model_*.cif")):
            k = int(cif.stem.rsplit("_", 1)[1])
            tag = f"DIAG_{label}_{t}_s{seed}_k{k}"
            mdl = score_prep.load(str(cif))
            mpdb, msdf = WD / f"{tag}_receptor.pdb", WD / f"{tag}_ligand.sdf"
            score_prep.write_pdb(mdl, mpdb, {"A"})
            err = score_prep.coords_to_sdf(score_prep.heavy(a for c in mdl if c.id == "L" for r in c for a in r), smi, msdf)
            conf = json.load(open(pdir / f"confidence_{t}_model_{k}.json"))
            row = dict(target=t, label=label, seed=seed, sample=k, job_rc=timing["rc"], wall_s=timing["wall_s"],
                       extra_flags=timing.get("extra_flags", ""), confidence_score=conf.get("confidence_score"),
                       ligand_iptm=conf.get("ligand_iptm"))
            if err:
                row["scoring_error"] = err
                rows.append(row)
                continue
            js, pb = OUT / f"{tag}.ostlig.json", OUT / f"{tag}.posebusters.csv"
            if not js.exists():
                score_tasks.run(score_tasks.OST + ["compare-ligand-structures", "-m", str(mpdb), "-ml", str(msdf),
                                                   "-r", str(rpdb), "-rl", str(rsdf), "--lddt-pli", "--rmsd", "-o", str(js)])
            if not pb.exists():
                score_tasks.run([score_tasks.BUST, str(msdf), "-p", str(mpdb), "--outfmt", "csv"], pb)
            score_parse.OUT = OUT
            m = score_parse.ligand_metrics(tag, 1)
            row.update(bisyrmsd=m.get("bisyrmsd"), lddt_pli=m.get("lddt_pli"), pb_valid=m.get("pb_valid"),
                       pb_failed_checks=";".join(k2.split(":", 1)[1] for k2 in m if k2.startswith("pb_fail:")),
                       success=m.get("success"))
            rows.append(row)
    # benchmark seed-1 Boltz-2 rows for the same targets, for side-by-side comparison
    sc = list(csv.DictReader(open(ROOT / "results/tables/scores.csv")))
    cf = list(csv.DictReader(open(ROOT / "results/tables/confidence.csv")))
    tm = {(r["target"], r["seed"]): r for r in csv.DictReader(open(ROOT / "results/tables/timing.csv"))
          if r["method"] == "boltz2"}
    for t in ["L1002", "L2001"]:
        for k in range(5):
            get = lambda met, src=sc: next((r["value"] for r in src if r["target"] == t and r["method"] == "boltz2"  # noqa: E731
                                            and r["seed"] == "1" and r["sample"] == str(k) and r["metric"] == met), None)
            fails = [r["metric"].split(":", 1)[1] for r in sc if r["target"] == t and r["method"] == "boltz2"
                     and r["seed"] == "1" and r["sample"] == str(k) and r["metric"].startswith("pb_fail:")]
            rows.append(dict(target=t, label="benchmark_potentials_off", seed=1, sample=k,
                             job_rc=tm[(t, "1")]["run_rc"], wall_s=tm[(t, "1")]["run_wall_s"], extra_flags="",
                             confidence_score=get("confidence_score", cf), ligand_iptm=get("ligand_iptm", cf),
                             bisyrmsd=get("bisyrmsd"), lddt_pli=get("lddt_pli"), pb_valid=get("pb_valid"),
                             pb_failed_checks=";".join(fails), success=get("success")))
    cols = ["target", "label", "seed", "sample", "bisyrmsd", "pb_valid", "pb_failed_checks", "success", "lddt_pli",
            "confidence_score", "ligand_iptm", "wall_s", "job_rc", "extra_flags", "scoring_error"]
    out = ROOT / "results/tables/diagnostics/boltz2_potentials.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["target"], r["label"], r["sample"])):
            w.writerow({c: r.get(c, "") for c in cols})
    print(f"{len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
