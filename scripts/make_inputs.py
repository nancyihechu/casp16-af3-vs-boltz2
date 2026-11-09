"""Build predictor inputs from data/targets/targets.json (and nothing else).

  python3 -I scripts/make_inputs.py msa     # inputs for the one-per-target MSA pass
  python3 -I scripts/make_inputs.py boltz   # Boltz YAMLs pointing at cached Boltz MSAs
  python3 -I scripts/make_inputs.py seeds   # per-seed AF3 JSONs from cached <name>_data.json

Rules enforced here:
  * Templates OFF. AF3: "templates": [] on every protein chain. Boltz: no `templates:`
    block. ColabFold: run without --templates (see run scripts).
  * Ligands are given as SMILES only. No native file is ever opened by this script.
  * Each method uses its OWN MSA pipeline. AF3: jackhmmer/nhmmer on the local AF3 DBs.
    Boltz-2 and ColabFold: the ColabFold MMseqs2 server, called through each tool's own
    code. MSAs are never shared between methods.

Seeds: the integers 1..5 for every method. NOTE: the same seed integer in different
codebases (JAX in AF3/ColabFold, PyTorch in Boltz) does NOT give comparable randomness.
Seeds are there for reproducibility only. The fairness control is the same NUMBER of
predictions per target per method: 5 seeds x 5 samples = 25 (AF3: 5 diffusion samples
per seed; Boltz: --diffusion_samples 5 per seed; ColabFold: 5 models per seed).
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TARGETS = json.load(open(ROOT / "data/targets/targets.json"))
INP = ROOT / "data/inputs"
MSA = ROOT / "msas"
SEEDS = [1, 2, 3, 4, 5]
AF2M_TARGETS = ["T1206", "H1215"]  # protein-only complexes


def af3_json(name, tgt, seeds):
    seqs = []
    for c in tgt["chains"]:
        ids = c["ids"] if len(c["ids"]) > 1 else c["ids"][0]
        if c["type"] == "protein":
            # templates [] + no MSA fields => MSA search runs, template search skipped
            seqs.append({"protein": {"id": ids, "sequence": c["sequence"], "templates": []}})
        elif c["type"] in ("rna", "dna"):
            seqs.append({c["type"]: {"id": ids, "sequence": c["sequence"]}})
        elif c["type"] == "ligand":
            seqs.append({"ligand": {"id": ids, "smiles": c["smiles"]}})
    return {"name": name, "modelSeeds": seeds, "dialect": "alphafold3", "version": 2, "sequences": seqs}


def boltz_protein_entities(tgt):
    """Unique protein sequences in input order, matching Boltz's entity grouping."""
    return [c for c in tgt["chains"] if c["type"] == "protein"]


def boltz_yaml(t, tgt):
    lines = ["version: 1", "sequences:"]
    for i, c in enumerate(tgt["chains"]):
        ids = "[" + ", ".join(c["ids"]) + "]"
        if c["type"] == "ligand":
            lines += ["  - ligand:", f"      id: {ids}", f"      smiles: '{c['smiles']}'"]
            continue
        lines += [f"  - {c['type']}:", f"      id: {ids}", f"      sequence: {c['sequence']}"]
        if c["type"] == "protein":
            pidx = boltz_protein_entities(tgt).index(c)
            msa = MSA / "boltz" / f"{t}_{pidx}.csv"
            assert msa.exists(), f"missing cached Boltz MSA {msa}; run the MSA pass first"
            lines.append(f"      msa: {msa}")
        # RNA/DNA: Boltz has no nucleic-acid MSA pipeline -> single sequence (its default)
    return "\n".join(lines) + "\n"


def colabfold_fasta(t, tgt):
    chains = [c["sequence"] for c in tgt["chains"] for _ in c["ids"]]
    return f">{t}\n{':'.join(chains)}\n"


def stage_msa():
    for d in ("af3_msa", "colabfold"):
        (INP / d).mkdir(parents=True, exist_ok=True)
    for t, tgt in TARGETS.items():
        json.dump(af3_json(t, tgt, [1]), open(INP / "af3_msa" / f"{t}.json", "w"), indent=2)
    for t in AF2M_TARGETS:
        (INP / "colabfold" / f"{t}.fasta").write_text(colabfold_fasta(t, TARGETS[t]))
    # Boltz MSA pass input: protein entities per target
    bz = {t: [c["sequence"] for c in boltz_protein_entities(tgt)] for t, tgt in TARGETS.items()}
    (INP / "boltz_msa_entities.json").write_text(json.dumps(bz, indent=2))
    print("wrote AF3 MSA-pass JSONs, ColabFold FASTAs, Boltz entity list")


def stage_seeds(only=None):
    (INP / "af3").mkdir(exist_ok=True)
    for t, tgt in TARGETS.items():
        if only and t not in only:
            continue
        data = MSA / "af3" / t.lower() / f"{t.lower()}_data.json"
        d = json.load(open(data))
        # sanity: the cached data.json must carry MSAs and zero templates for every protein chain
        for s in d["sequences"]:
            if "protein" in s:
                p = s["protein"]
                assert p["templates"] == [], f"{t}: templates present in cached data.json"
                assert p["unpairedMsa"] is not None and p["pairedMsa"] is not None
        for seed in SEEDS:
            d2 = dict(d, name=f"{t}_s{seed}", modelSeeds=[seed])
            json.dump(d2, open(INP / "af3" / f"{t}_s{seed}.json", "w"))
    print("wrote per-seed AF3 JSONs (from cached data.json) for", only or "all targets")


def stage_boltz():
    (INP / "boltz").mkdir(exist_ok=True)
    for t, tgt in TARGETS.items():
        (INP / "boltz" / f"{t}.yaml").write_text(boltz_yaml(t, tgt))
    print("wrote Boltz YAMLs (with cached MSAs)")


if __name__ == "__main__":
    fn = {"msa": stage_msa, "seeds": stage_seeds, "boltz": stage_boltz}[sys.argv[1]]
    fn(sys.argv[2:]) if sys.argv[1] == "seeds" and len(sys.argv) > 2 else fn()
