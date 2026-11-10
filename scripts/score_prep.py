"""Stage D prep: build the scoring task list and the ligand SDF / receptor PDB files.

  envs/p3_tools/bin/python -I scripts/score_prep.py

Reads results/tables/manifest.csv (run collect_runs.py first) and the natives in
data/targets/natives/. Natives are used ONLY here and in scoring, never by a predictor.

Ligands (L1002, L2001): OpenStructure needs ligand connectivity, so every ligand
(reference and model) is written to SDF. Connectivity comes from coordinates
(rdDetermineBonds.DetermineConnectivity), and bond orders from the official CASP
SMILES (AssignBondOrdersFromTemplate). A ligand whose conversion fails is recorded
as a scoring failure, never silently dropped.

M1211 RNA-only pass: model chain C and native chain 0 are written out alone and scored
with the same OST TM-score / lDDT as R1211 (task kind "rna_only").

L2001 reference: CASP stoichiometry is A1, but the crystal ASU has 2 protein chains.
We score against ONE chain: the protein chain closest to the reference ligand.

Output: scoring/tasks.tsv plus scoring/work/** (git-ignored with predictions/).
"""
import csv
import json
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds

ROOT = Path(__file__).resolve().parent.parent
NAT = ROOT / "data/targets/natives"
WORK = ROOT / "predictions/scoring_work"
TARGETS = json.load(open(ROOT / "data/targets/targets.json"))
LIG = {"L1002": ("L1002_protein.pdb", "L1002_ligand_293.pdb"),
       "L2001": ("L2001_protein.pdb", "L2001_ligand_761.pdb")}
NATIVE = {t: f"{t}.pdb" for t in ["T1206", "H1215", "M1211", "R1211", "M1276"]}


def load(path):
    p = MMCIFParser(QUIET=True) if str(path).endswith(".cif") else PDBParser(QUIET=True)
    return next(iter(p.get_structure("s", str(path))))


def coords_to_sdf(atoms, smiles, out):
    """atoms: list of (element, xyz). Returns None or an error string."""
    try:
        rw = Chem.RWMol()
        conf = Chem.Conformer(len(atoms))
        for i, (el, xyz) in enumerate(atoms):
            rw.AddAtom(Chem.Atom(el.capitalize()))
            conf.SetAtomPosition(i, [float(v) for v in xyz])
        m = rw.GetMol()
        m.AddConformer(conf, assignId=True)
        rdDetermineBonds.DetermineConnectivity(m)
        tmpl = Chem.MolFromSmiles(smiles)
        m = AllChem.AssignBondOrdersFromTemplate(tmpl, m)
        # Atoms built from bare coordinates carry radical flags / no implicit H after the
        # template step. Clear them so PoseBusters does not flag a conversion artefact as a
        # model error (hydrogen count is then implied by valence, as in the SMILES).
        for a in m.GetAtoms():
            a.SetNumRadicalElectrons(0)
            a.SetNoImplicit(False)
        Chem.SanitizeMol(m)
        assert Chem.MolToInchiKey(Chem.RemoveHs(m)).split("-")[0] == Chem.MolToInchiKey(tmpl).split("-")[0], \
            "connectivity/formula does not match the SMILES"
        Chem.AssignStereochemistryFrom3D(m)
        Chem.MolToMolFile(m, str(out))
        return None
    except Exception as e:  # noqa: BLE001
        return f"{type(e).__name__}: {e}"


def heavy(atoms):
    return [(a.element, a.coord) for a in atoms if a.element != "H"]


class ChainOnly(Select):
    """Polymer residues of one chain (used for the M1211 RNA-only pass)."""
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, c):
        return c.id == self.chain_id

    def accept_residue(self, r):
        return r.id[0] == " "


class ProteinOnly(Select):
    def __init__(self, chain_ids=None):
        self.chain_ids = chain_ids

    def accept_chain(self, c):
        return self.chain_ids is None or c.id in self.chain_ids

    def accept_residue(self, r):
        return r.id[0] == " "  # standard polymer residues only


def write_pdb(model, out, chain_ids=None):
    io = PDBIO()
    io.set_structure(model)
    io.save(str(out), ProteinOnly(chain_ids))


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    tasks, failures = [], []
    # Reference ligand targets: SDF + single receptor chain
    ref = {}
    for t, (prot, lig) in LIG.items():
        smi = next(c["smiles"] for c in TARGETS[t]["chains"] if c["type"] == "ligand")
        lig_atoms = heavy(load(NAT / lig).get_atoms())
        rsdf = WORK / f"ref_{t}_ligand.sdf"
        err = coords_to_sdf(lig_atoms, smi, rsdf)
        assert err is None, f"reference ligand {t} failed: {err}"
        pm = load(NAT / prot)
        lc = np.array([xyz for _, xyz in lig_atoms])
        dist = {c.id: min(np.linalg.norm(np.array([a.coord for a in c.get_atoms()])[:, None] - lc[None], axis=2).min(), 1e9)
                for c in pm}
        chain = min(dist, key=dist.get)
        rpdb = WORK / f"ref_{t}_receptor_chain{chain}.pdb"
        write_pdb(pm, rpdb, {chain})
        ref[t] = (rsdf, rpdb, smi)
        print(f"{t}: reference receptor chain {chain} (min dist to ligand {dist})")

    rna_ref = WORK / "ref_M1211_rna_chain0.pdb"
    for r in csv.DictReader(open(ROOT / "results/tables/manifest.csv")):
        t, m, s, k, path = r["target"], r["method"], r["seed"], r["sample"], r["path"]
        tag = f"{t}_{m}_s{s}_k{k}"
        if t in LIG:
            rsdf, rpdb, smi = ref[t]
            mdl = load(path)
            lig_atoms = heavy(a for c in mdl for res in c if res.id[0] != " " or c.id == "L" for a in res)
            msdf, mpdb = WORK / f"{tag}_ligand.sdf", WORK / f"{tag}_receptor.pdb"
            write_pdb(mdl, mpdb, {"A"})
            err = coords_to_sdf(lig_atoms, smi, msdf)
            if err:
                failures.append((tag, err))
                tasks.append([tag, t, m, s, k, "ligand_FAILED_SDF", path, "", "", "", ""])
                continue
            tasks.append([tag, t, m, s, k, "ligand", str(mpdb), str(msdf), str(rpdb), str(rsdf), path])
        else:
            tasks.append([tag, t, m, s, k, "polymer", path, "", str(NAT / NATIVE[t]), "", ""])
            if t == "M1211":
                # RNA chain alone (model chain C vs native RNA chain 0), scored like R1211
                if not rna_ref.exists():
                    io = PDBIO(); io.set_structure(load(NAT / "M1211.pdb")); io.save(str(rna_ref), ChainOnly("0"))
                rp = WORK / f"{tag}_rnaonly.pdb"
                io = PDBIO(); io.set_structure(load(path)); io.save(str(rp), ChainOnly("C"))
                tasks.append([f"{tag}_rnaonly", t, m, s, k, "rna_only", str(rp), "", str(rna_ref), "", path])

    with open(ROOT / "predictions/scoring_work/tasks.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["tag", "target", "method", "seed", "sample", "kind", "model", "model_ligand",
                    "reference", "reference_ligand", "model_full"])
        w.writerows(tasks)
    print(f"{len(tasks)} tasks; {len(failures)} ligand SDF failures")
    for f in failures:
        print("  FAIL", *f)


if __name__ == "__main__":
    main()
