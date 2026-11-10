"""Walk predictions/ and write one row per prediction (the run manifest).

  envs/p3_tools/bin/python -I scripts/collect_runs.py [--include-positive-control]

Output: results/tables/manifest.csv with columns
  target, method, seed, sample, path, run_rc, run_wall_s, run_inference_s, msa_s, msa_source,
  chains (e.g. "A:protein:237;B:protein:237;L:ligand:1"), and confidence columns
  (ptm, iptm, ranking_score, confidence_score, plddt_mean, ...; blank where a method has none).

MSA time is per target per method (one cached MSA pass) and comes from msas/msa_timing_*.tsv.
It is not split across seeds; run_wall_s covers the per-seed run with the MSA already cached.
"""
import csv
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

ROOT = Path(__file__).resolve().parent.parent
PRED = ROOT / "predictions"
NUC = {"A", "C", "G", "U", "DA", "DC", "DG", "DT"}
AA3 = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())


def chains_of(path):
    p = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    s = p.get_structure("m", path)
    out = []
    for ch in next(iter(s)):
        res = list(ch)
        names = {r.get_resname().strip() for r in res}
        if names & AA3:
            kind = "protein"
        elif names <= NUC:
            kind = "dna" if any(n.startswith("D") for n in names) else "rna"
        else:
            kind = "ligand"
        out.append(f"{ch.id}:{kind}:{len(res)}")
    return ";".join(out)


def read_timing(p):
    if not Path(p).exists():
        return {}
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    return rows[0] if rows else {}


def msa_times():
    m = {}
    for row in csv.DictReader(open(ROOT / "msas/msa_timing_af3.tsv"), delimiter="\t"):
        m[(row["target"], "af3")] = (row["data_pipeline_s"], "AF3 jackhmmer/nhmmer, local AF3 DBs")
    for row in csv.DictReader(open(ROOT / "msas/msa_timing_boltz.tsv"), delimiter="\t"):
        m[(row["target"], "boltz2")] = (row["msa_s"], "ColabFold MMseqs2 server via boltz compute_msa")
    for row in csv.DictReader(open(ROOT / "msas/msa_timing_colabfold.tsv"), delimiter="\t"):
        m[(row["target"], "af2m")] = (row["msa_s"], "ColabFold MMseqs2 server via colabfold_batch")
    return m


def af3_rows():
    for d in sorted(PRED.glob("af3/*/s*")):
        t, seed = d.parent.name, int(d.name[1:])
        name = f"{t}_s{seed}".lower()
        tm = read_timing(d / f"{t}_s{seed}_timing.tsv")
        for sd in sorted((d / name).glob("seed-*_sample-*")):
            k = int(sd.name.split("sample-")[1])
            sc = json.load(open(sd / "summary_confidences.json"))
            conf = json.load(open(sd / "confidences.json"))
            yield dict(target=t, method="af3", seed=seed, sample=k, path=str(sd / "model.cif"),
                       run_rc=tm.get("rc"), run_wall_s=tm.get("wall_s"), run_inference_s=tm.get("inference_s"),
                       chains=chains_of(str(sd / "model.cif")),
                       ptm=sc.get("ptm"), iptm=sc.get("iptm"), ranking_score=sc.get("ranking_score"),
                       plddt_mean=round(float(np.mean(conf["atom_plddts"])), 2),
                       has_clash=sc.get("has_clash"), fraction_disordered=sc.get("fraction_disordered"),
                       chain_pair_iptm=json.dumps(sc.get("chain_pair_iptm")))


def boltz_rows():
    for d in sorted(PRED.glob("boltz2/*/s*")):
        t, seed = d.parent.name, int(d.name[1:])
        tm = read_timing(d / "timing.tsv")
        pdir = d / f"boltz_results_{t}" / "predictions" / t
        for cif in sorted(pdir.glob(f"{t}_model_*.cif")):
            k = int(cif.stem.rsplit("_", 1)[1])
            c = json.load(open(pdir / f"confidence_{t}_model_{k}.json"))
            yield dict(target=t, method="boltz2", seed=seed, sample=k, path=str(cif),
                       run_rc=tm.get("rc"), run_wall_s=tm.get("wall_s"), run_inference_s=tm.get("inference_s"),
                       chains=chains_of(str(cif)),
                       ptm=c.get("ptm"), iptm=c.get("iptm"), confidence_score=c.get("confidence_score"),
                       plddt_mean=round(100 * c["complex_plddt"], 2) if "complex_plddt" in c else None,
                       ligand_iptm=c.get("ligand_iptm"), protein_iptm=c.get("protein_iptm"),
                       complex_iplddt=c.get("complex_iplddt"),
                       chain_pair_iptm=json.dumps(c.get("pair_chains_iptm")))


def af2m_seed_inference(logfile):
    """Per-seed inference = sum of ColabFold's own per-model 'took Xs' lines for that seed."""
    out = {}
    for line in open(logfile):
        m = re.search(r"model_\d_seed_(\d+) took ([0-9.]+)s", line)
        if m:
            out[int(m.group(1))] = out.get(int(m.group(1)), 0.0) + float(m.group(2))
    return out


def af2m_rows():
    for d in sorted(PRED.glob("af2m/*/s*_n*")):
        t = d.parent.name
        tm = read_timing(d / "timing.tsv")
        n_in_job = int(d.name.split("_n")[1])
        inf = af2m_seed_inference(d / "colabfold.log")
        for js in sorted(d.glob(f"{t}_scores_rank_*_model_*_seed_*.json")):
            m = re.search(r"model_(\d)_seed_(\d+)", js.name)
            model, seed = int(m.group(1)), int(m.group(2))
            pdb = Path(str(js).replace("_scores_", "_unrelaxed_").replace(".json", ".pdb"))
            sc = json.load(open(js))
            yield dict(target=t, method="af2m", seed=seed, sample=model - 1, path=str(pdb),
                       run_rc=tm.get("rc"), run_wall_s=tm.get("wall_s"), run_inference_s=round(inf[seed], 1),
                       seeds_in_job=n_in_job,
                       chains=chains_of(str(pdb)), ptm=sc.get("ptm"), iptm=sc.get("iptm"),
                       ranking_score=round(0.8 * sc["iptm"] + 0.2 * sc["ptm"], 4),
                       plddt_mean=round(float(np.mean(sc["plddt"])), 2))


def main():
    mt = msa_times()
    rows = [*af3_rows(), *boltz_rows(), *af2m_rows()]
    for r in rows:
        r["msa_s"], r["msa_source"] = mt.get((r["target"], r["method"]), ("", ""))
        r.setdefault("seeds_in_job", 1)
    cols = ["target", "method", "seed", "sample", "path", "run_rc", "run_wall_s", "seeds_in_job", "run_inference_s",
            "msa_s", "msa_source", "chains", "ptm", "iptm", "ranking_score", "confidence_score",
            "plddt_mean", "complex_iplddt", "ligand_iptm", "protein_iptm", "has_clash",
            "fraction_disordered", "chain_pair_iptm"]
    out = ROOT / "results/tables/manifest.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["target"], r["method"], r["seed"], r["sample"])):
            w.writerow({k: r.get(k) for k in cols})
    print(f"{len(rows)} predictions -> {out}")


if __name__ == "__main__":
    main()
