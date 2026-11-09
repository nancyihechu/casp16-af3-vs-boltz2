"""Boltz-2 MSA pass: one call per target to Boltz's OWN `boltz.main.compute_msa`.

This is the same function `boltz predict --use_msa_server` calls internally, with the same
defaults (ColabFold MMseqs2 server https://api.colabfold.com, pairing strategy 'greedy',
use_env=True). Pairing happens only when a target has >1 protein entity, exactly as in
`boltz predict`. Output: msas/boltz/<T>_<i>.csv, referenced later via the YAML `msa:` field,
so `boltz predict` skips MSA generation for every seed.

Boltz builds no MSA for RNA/DNA/ligands (it sets msa_id=-1 for non-protein chains).

Run with envs/boltz2/bin/python -I scripts/boltz_msa.py  (needs outbound network).
"""
import json
import time
from pathlib import Path

from boltz.main import compute_msa

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "msas" / "boltz"
OUT.mkdir(parents=True, exist_ok=True)
ents = json.load(open(ROOT / "data/inputs/boltz_msa_entities.json"))
URL, PAIRING = "https://api.colabfold.com", "greedy"  # boltz predict defaults

with open(ROOT / "msas" / "msa_timing_boltz.tsv", "w") as tsv:
    tsv.write("target\tmethod\tn_protein_entities\tmsa_s\tdepths\n")
    for t, seqs in ents.items():
        if not seqs:
            tsv.write(f"{t}\tboltz2\t0\t0\tno protein entity (no MSA in Boltz)\n")
            print(t, "no protein entities; skipped")
            continue
        data = {f"{t}_{i}": s for i, s in enumerate(seqs)}
        t0 = time.time()
        compute_msa(data=data, target_id=t, msa_dir=OUT, msa_server_url=URL, msa_pairing_strategy=PAIRING)
        dt = time.time() - t0
        depths = []
        for k in data:
            rows = (OUT / f"{k}.csv").read_text().splitlines()[1:]
            npair = sum(1 for r in rows if not r.startswith("-1,"))
            depths.append(f"{k}:{len(rows)}(paired_rows={npair})")
        tsv.write(f"{t}\tboltz2\t{len(seqs)}\t{dt:.1f}\t{';'.join(depths)}\n")
        tsv.flush()
        print(t, f"{dt:.1f}s", depths)
