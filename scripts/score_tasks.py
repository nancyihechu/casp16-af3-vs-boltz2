"""Stage D: score every prediction in predictions/scoring_work/tasks.tsv.

  envs/p3_tools/bin/python -I scripts/score_tasks.py <n_workers> [<shard> <n_shards>]

Per task (raw outputs go to predictions/scoring_work/out/<tag>.*):
  polymer: OpenStructure 2.12.0 `compare-structures` with --lddt --ilddt --dockq --ics --ips
           --qs-score --rigid-scores --tm-score (built-in USalign). Chain mapping and residue
           alignment are done by OST (common residues only).
           Protein-protein targets also get standalone `USalign -mm 1 -ter 0` (as requested).
  rna_only (M1211 chain C vs native RNA chain 0): OST compare-structures --lddt --tm-score,
           i.e. the same scorer and options that give R1211 its tm_score / lddt.
  ligand:  OpenStructure `compare-ligand-structures --lddt-pli --rmsd` (BiSyRMSD, the CASP16
           ligand metric). PoseBusters 0.6.5 in `dock` mode (predicted ligand + the model's own
           protein). Its redock-mode `rmsd_<=_2A` check is unusable here because model and
           reference are in different frames, so RMSD comes only from OST BiSyRMSD.
A failed scorer writes <tag>.err. It is never retried with different settings.
"""
import csv
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
W = ROOT / "predictions/scoring_work"
OUT = W / "out"
OUT.mkdir(parents=True, exist_ok=True)
NANCY = ROOT.parent
OST = ["singularity", "exec", "-B", "/cluster:/cluster", str(NANCY / "containers/ost.sif"), "ost"]
USALIGN = str(NANCY / "envs/p3_tools/bin/USalign")
BUST = str(NANCY / "envs/p3_tools/bin/bust")
PP = {"T1206", "H1215"}


def run(cmd, out_txt=None):
    p = subprocess.run(cmd, capture_output=True, text=True)
    if out_txt:
        Path(out_txt).write_text(p.stdout)
    return p.returncode, p.stderr[-3000:]


def score(t):
    tag, errs = t["tag"], []
    if t["kind"] == "polymer":
        js = OUT / f"{tag}.ost.json"
        if not js.exists():
            rc, e = run(OST + ["compare-structures", "-m", t["model"], "-r", t["reference"], "-o", str(js),
                               "--lddt", "--ilddt", "--dockq", "--ics", "--ips", "--qs-score",
                               "--rigid-scores", "--tm-score"])
            if rc:
                errs.append(f"ost rc={rc}: {e}")
        if t["target"] in PP:
            us = OUT / f"{tag}.usalign.txt"
            if not us.exists():
                rc, e = run([USALIGN, t["model"], t["reference"], "-mm", "1", "-ter", "0", "-outfmt", "2"], us)
                if rc:
                    errs.append(f"usalign rc={rc}: {e}")
    elif t["kind"] == "rna_only":
        js = OUT / f"{tag}.ost.json"
        if not js.exists():
            rc, e = run(OST + ["compare-structures", "-m", t["model"], "-r", t["reference"], "-o", str(js),
                               "--lddt", "--tm-score"])
            if rc:
                errs.append(f"ost rc={rc}: {e}")
    elif t["kind"] == "ligand":
        js = OUT / f"{tag}.ostlig.json"
        if not js.exists():
            rc, e = run(OST + ["compare-ligand-structures", "-m", t["model"], "-ml", t["model_ligand"],
                               "-r", t["reference"], "-rl", t["reference_ligand"], "--lddt-pli", "--rmsd",
                               "-o", str(js)])
            if rc:
                errs.append(f"ost-ligand rc={rc}: {e}")
        pb = OUT / f"{tag}.posebusters.csv"
        if not pb.exists():
            rc, e = run([BUST, t["model_ligand"], "-p", t["model"], "--outfmt", "csv"], pb)
            if rc:
                errs.append(f"posebusters rc={rc}: {e}")
    else:
        errs.append(f"not scorable: {t['kind']}")
    if errs:
        (OUT / f"{tag}.err").write_text("\n".join(errs))
    return tag, not errs


def main():
    nw = int(sys.argv[1])
    shard, nsh = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, 1)
    tasks = list(csv.DictReader(open(W / "tasks.tsv"), delimiter="\t"))[shard::nsh]
    ok = 0
    with ThreadPoolExecutor(nw) as ex:
        for tag, good in ex.map(score, tasks):
            ok += good
            if not good:
                print("FAILED", tag, flush=True)
    print(f"shard {shard}/{nsh}: {ok}/{len(tasks)} scored OK")


if __name__ == "__main__":
    main()
