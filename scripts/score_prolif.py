"""Stage D (part 3c): protein-ligand interaction recovery with ProLIF.

  envs/p3_tools/bin/python -I scripts/score_prolif.py <n_workers> [tag ...]

For every ligand task in predictions/scoring_work/tasks.tsv (run score_prep.py and
score_tasks.py first; the OST ligand output is needed for the model->reference ligand
assignment on multi-ligand targets):

  1. Hydrogens, IDENTICAL procedure for native and model:
     - protein: all existing H stripped (CASP natives ship with H), then PDBFixer 1.12.0 (OpenMM 8.3.1) addMissingHydrogens(pH=7.0). Missing residues
       and missing heavy atoms are NOT added (cleared before protonation), so only H change.
     - ligand: the SMILES-bond-ordered SDF from score_prep.py, RDKit Chem.AddHs(addCoords=True).
     - protein bonds = OpenMM residue-template bonds written as explicit CONECT records, read
       by RDKit with proximityBonding=False (no distance-based bonds).
  2. Fingerprint: ProLIF 2.2.3, plf.Fingerprint() with its DEFAULT interaction set, count=False.
     An interaction = (protein residue, interaction type).
  3. Residue identity: model and native numbering differ, so every model residue is mapped to
     a native residue by global sequence alignment within each chain pair (Biopython
     PairwiseAligner). Chain pairs: identity for single-chain targets; for L4001 the chain
     mapping OST used for that ligand.
  4. recovery = |native interactions reproduced by the pose| / |native interactions|.
     Multi-ligand targets (L4001): each model ligand is compared with the reference ligand OST
     assigned it to; native interactions of an unassigned reference ligand count as not
     recovered; the target value pools both copies (sum recovered / sum native).

This is NOT PoseBench PLIF-WM (see CLAUDE.md): no histogram / Wasserstein step and no
dataset-relative normalisation.
Output: predictions/scoring_work/out/<tag>.prolif.json
"""
import csv
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

warnings.filterwarnings("ignore")
import prolif as plf  # noqa: E402
from Bio.Align import PairwiseAligner  # noqa: E402
from Bio.PDB import PDBParser  # noqa: E402
from Bio.SeqUtils import seq1  # noqa: E402
from openmm.app import PDBFile  # noqa: E402
from pdbfixer import PDBFixer  # noqa: E402
from rdkit import Chem  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
W = ROOT / "predictions/scoring_work"
OUT = W / "out"
HDIR = W / "protonated"
HDIR.mkdir(exist_ok=True)


def strip_h(pdb_in):
    """Remove every hydrogen first, so native (CASP files already carry H) and model end up
    with H from the same tool."""
    noh = HDIR / (Path(pdb_in).stem + "_noH.pdb")
    keep = [l for l in open(pdb_in) if not (l.startswith(("ATOM", "HETATM")) and
                                             (l[76:78].strip() == "H" or (l[76:78].strip() == "" and l[12:16].strip().startswith("H"))))]
    noh.write_text("".join(keep))
    return noh


def protonate(pdb_in):
    out = HDIR / (Path(pdb_in).stem + "_H.pdb")
    if not out.exists():
        fx = PDBFixer(filename=str(strip_h(pdb_in)))
        fx.findMissingResidues()
        fx.missingResidues = {}
        fx.findMissingAtoms()
        fx.missingAtoms, fx.missingTerminals = {}, {}
        fx.addMissingHydrogens(7.0)
        tmp = HDIR / (Path(pdb_in).stem + "_H_noconect.pdb")
        with open(tmp, "w") as f:
            PDBFile.writeFile(fx.topology, fx.positions, f, keepIds=True)
        # Explicit CONECT for EVERY bond in OpenMM's template-based topology, so RDKit reads
        # bonds from residue templates (not from distances). Written ATOM/HETATM records follow
        # topology.atoms() order, which gives the atom-index -> serial map.
        lines = [l for l in open(tmp) if not l.startswith(("CONECT", "END"))]
        serials = [int(l[6:11]) for l in lines if l.startswith(("ATOM", "HETATM"))]
        atoms = list(fx.topology.atoms())
        assert len(serials) == len(atoms), "atom count mismatch writing CONECT"
        nb = {}
        for a, b in fx.topology.bonds():
            nb.setdefault(a.index, []).append(serials[b.index])
            nb.setdefault(b.index, []).append(serials[a.index])
        con = []
        for i, partners in sorted(nb.items()):
            for k in range(0, len(partners), 4):
                con.append("CONECT" + f"{serials[i]:5d}" + "".join(f"{x:5d}" for x in partners[k:k + 4]) + "\n")
        out.write_text("".join(lines) + "".join(con) + "END\n")
        tmp.unlink()
    return out


def prot_mol(pdb_h):
    # proximityBonding=False: bonds come only from the CONECT records written in protonate(),
    # i.e. OpenMM residue templates (+ disulfides), so crystal close contacts (e.g. L2001 Arg170
    # NH2 - Arg174 NH1, 1.85 A) do not become bonds.
    m = Chem.MolFromPDBFile(str(pdb_h), removeHs=False, proximityBonding=False)
    if m is None:
        raise ValueError(f"RDKit could not read {pdb_h}")
    return plf.Molecule(m)


def lig_mol(sdf):
    m = Chem.MolFromMolFile(str(sdf))
    return plf.Molecule.from_rdkit(Chem.AddHs(m, addCoords=True))


def interactions(lig_sdf, prot):
    fp = plf.Fingerprint(count=False)
    fp.run_from_iterable([lig_mol(lig_sdf)], prot, progress=False, n_jobs=1)
    ifp = fp.ifp[0]
    out = set()
    for (_lres, pres), inter in ifp.items():
        for itype in inter:
            out.add((pres.chain, pres.number, pres.name, itype))
    return out, sorted(fp.interactions)


def chain_seqs(pdb):
    s = next(iter(PDBParser(QUIET=True).get_structure("s", str(pdb))))
    return {c.id: [(r.id[1], seq1(r.get_resname())) for r in c if r.id[0] == " "] for c in s}


def residue_map(mdl_pdb, ref_pdb, chain_pairs):
    """{(mdl_chain, mdl_resnum): (ref_chain, ref_resnum)} via global alignment per chain pair."""
    al = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                         open_gap_score=-5, extend_gap_score=-0.5)
    ms, rs = chain_seqs(mdl_pdb), chain_seqs(ref_pdb)
    mp = {}
    for mc, rc in chain_pairs.items():
        a = al.align("".join(x[1] for x in ms[mc]), "".join(x[1] for x in rs[rc]))[0]
        for (m0, m1), (r0, r1) in zip(*a.aligned):
            for i, j in zip(range(m0, m1), range(r0, r1)):
                mp[(mc, ms[mc][i][0])] = (rc, rs[rc][j][0])
    return mp


def score(t):
    tag = t["tag"]
    out = OUT / f"{tag}.prolif.json"
    if out.exists():
        return tag, True, ""
    try:
        ref_ligs = t["reference_ligand"].split(";")
        mdl_ligs = t["model_ligand"].split(";")
        ref_h, mdl_h = protonate(t["reference"]), protonate(t["model"])
        ref_p, mdl_p = prot_mol(ref_h), prot_mol(mdl_h)
        # model ligand -> reference ligand and chain mapping, from OST
        if len(ref_ligs) == 1:
            ref_chain = next(iter(chain_seqs(t["reference"])))
            pairs = [(mdl_ligs[0], ref_ligs[0], {next(iter(chain_seqs(t["model"]))): ref_chain})]
        else:
            ost = json.load(open(OUT / f"{tag}.ostlig.json"))
            pairs = []
            for sc in ost["rmsd"]["assigned_scores"]:
                cm = sc["chain_mapping"]  # OST: {reference_chain: model_chain}
                pairs.append((sc["model_ligand"], sc["reference_ligand"], {v: k for k, v in cm.items()}))
        res = {"per_ligand": [], "hydrogens": "protein PDBFixer addMissingHydrogens pH7.0 (no heavy atoms added); "
                                                "ligand RDKit AddHs(addCoords=True) on SMILES-bonded SDF"}
        tot_nat = tot_rec = 0
        assigned_refs = set()
        for mlig, rlig, cpairs in pairs:
            nat, itypes = interactions(rlig, ref_p)
            mdl_raw, _ = interactions(mlig, mdl_p)
            rmap = residue_map(t["model"], t["reference"], cpairs)
            mdl = set()
            for (c, n, name, it) in mdl_raw:
                if (c, n) in rmap:
                    rc, rn = rmap[(c, n)]
                    mdl.add((rc, rn, it))
            nat_k = {(c, n, it) for (c, n, _nm, it) in nat}
            rec = len(nat_k & mdl)
            res["per_ligand"].append({"model_ligand": mlig, "reference_ligand": rlig, "n_native": len(nat_k),
                                      "n_model": len(mdl_raw), "n_recovered": rec,
                                      "recovery": rec / len(nat_k) if nat_k else None})
            tot_nat += len(nat_k)
            tot_rec += rec
            assigned_refs.add(rlig)
            res["interaction_types"] = itypes
        for rlig in ref_ligs:  # unassigned reference ligands: all their interactions are missed
            if rlig not in assigned_refs:
                nat, _ = interactions(rlig, ref_p)
                n = len({(c, n, it) for (c, n, _nm, it) in nat})
                res["per_ligand"].append({"model_ligand": None, "reference_ligand": rlig, "n_native": n,
                                          "n_model": 0, "n_recovered": 0, "recovery": 0.0})
                tot_nat += n
        res["n_native"], res["n_recovered"] = tot_nat, tot_rec
        res["recovery"] = tot_rec / tot_nat if tot_nat else None
        out.write_text(json.dumps(res, indent=1))
        return tag, True, ""
    except Exception as e:  # noqa: BLE001
        (OUT / f"{tag}.prolif.err").write_text(f"{type(e).__name__}: {e}")
        return tag, False, f"{type(e).__name__}: {e}"


def main():
    nw = int(sys.argv[1])
    tasks = [t for t in csv.DictReader(open(W / "tasks.tsv"), delimiter="\t") if t["kind"] == "ligand"]
    if len(sys.argv) > 2:
        tasks = [t for t in tasks if t["tag"] in sys.argv[2:]]
    for ref in sorted({t["reference"] for t in tasks}):  # shared files: protonate once, serially
        protonate(ref)
    ok = 0
    with ProcessPoolExecutor(nw) as ex:
        for tag, good, err in ex.map(score, tasks):
            ok += good
            if not good:
                print("FAILED", tag, err, flush=True)
    print(f"ProLIF: {ok}/{len(tasks)} ligand predictions scored")


if __name__ == "__main__":
    main()
