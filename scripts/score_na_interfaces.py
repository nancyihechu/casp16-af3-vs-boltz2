"""Stage D (part 3b): interface lDDT per protein-nucleic-acid interface (M1211, M1276).

  envs/p3_tools/bin/python -I scripts/score_na_interfaces.py <n_workers>

OST 2.12.0 reports `ilddt` only over ALL inter-chain contacts (for M1211 that includes the
protein-protein A-B interface). To restrict to protein-NA interfaces, the same way
dockq_protein_na_mean is restricted, each protein-NA reference interface (r1, r2) listed in
`contact_reference_interfaces` of the full-complex OST run is re-scored on its own:
  - reference: only chains r1, r2 of the native
  - model: only the model chains OST mapped to them in the full-complex run (chain_mapping is
    {reference_chain: model_chain})
  - `ost compare-structures --ilddt -c r1:m1 r2:m2` (mapping enforced, identical to the full run)
Output: predictions/scoring_work/out/<tag>.naiface_<r1>-<r2>.json
ICS per interface needs no re-scoring: it is `per_interface_ics` in the full-complex output.
"""
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select

ROOT = Path(__file__).resolve().parent.parent
W = ROOT / "predictions/scoring_work"
OUT = W / "out"
SUB = W / "naiface"
SUB.mkdir(exist_ok=True)
OST = ["singularity", "exec", "-B", "/cluster:/cluster", str(ROOT.parent / "containers/ost.sif"), "ost"]
NA_REF_CHAINS = {"M1211": {"0"}, "M1276": {"0"}}


class Chains(Select):
    def __init__(self, ids):
        self.ids = set(ids)

    def accept_chain(self, c):
        return c.id in self.ids

    def accept_residue(self, r):
        return r.id[0] == " "


def subset(path, ids, out):
    p = MMCIFParser(QUIET=True) if path.endswith(".cif") else PDBParser(QUIET=True)
    io = PDBIO()
    io.set_structure(next(iter(p.get_structure("s", path))))
    io.save(str(out), Chains(ids))


def work(t):
    tag, tg = t["tag"], t["target"]
    full = json.load(open(OUT / f"{tag}.ost.json"))
    cm = full["chain_mapping"]
    res = []
    for r1, r2 in full["contact_reference_interfaces"]:
        if (r1 in NA_REF_CHAINS[tg]) == (r2 in NA_REF_CHAINS[tg]):
            continue  # not a protein-NA interface
        js = OUT / f"{tag}.naiface_{r1}-{r2}.json"
        if js.exists():
            res.append(True)
            continue
        if r1 not in cm or r2 not in cm:
            js.write_text(json.dumps({"status": "UNMAPPED", "chain_mapping": cm}))
            res.append(False)
            continue
        m1, m2 = cm[r1], cm[r2]
        rp, mp = SUB / f"{tg}_ref_{r1}{r2}.pdb", SUB / f"{tag}_{m1}{m2}.pdb"
        if not rp.exists():
            subset(t["reference"], [r1, r2], rp)
        subset(t["model"], [m1, m2], mp)
        p = subprocess.run(OST + ["compare-structures", "-m", str(mp), "-r", str(rp), "-o", str(js),
                                  "--ilddt", "-c", f"{r1}:{m1}", f"{r2}:{m2}"], capture_output=True, text=True)
        if p.returncode:
            (OUT / f"{tag}.naiface_{r1}-{r2}.err").write_text(p.stderr[-3000:])
        res.append(p.returncode == 0)
    return tag, all(res), len(res)


def main():
    nw = int(sys.argv[1])
    tasks = [t for t in csv.DictReader(open(W / "tasks.tsv"), delimiter="\t")
             if t["kind"] == "polymer" and t["target"] in NA_REF_CHAINS]
    # reference subsets are shared: build them serially first
    for t in tasks:
        full = json.load(open(OUT / f"{t['tag']}.ost.json"))
        for r1, r2 in full["contact_reference_interfaces"]:
            rp = SUB / f"{t['target']}_ref_{r1}{r2}.pdb"
            if (r1 in NA_REF_CHAINS[t["target"]]) != (r2 in NA_REF_CHAINS[t["target"]]) and not rp.exists():
                subset(t["reference"], [r1, r2], rp)
    ok = n = 0
    with ThreadPoolExecutor(nw) as ex:
        for tag, good, k in ex.map(work, tasks):
            ok += good
            n += k
            if not good:
                print("FAILED", tag, flush=True)
    print(f"{ok}/{len(tasks)} predictions; {n} protein-NA interfaces scored")


if __name__ == "__main__":
    main()
