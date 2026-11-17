"""Build data/targets/targets.json: the single source of truth for every predictor input.

Sequences come ONLY from official CASP16 files in data/targets/casp16_downloads/
(sequences/ folder, or the L1000/L2000 target pages for the ligand proteins).
SMILES come ONLY from the official CASP SMILES .tsv files.
Native structures / ligand_*.pdb are never read here.

Stoichiometry follows CASP (see CLAUDE.md), not the chain count of the native file.

Run: python3 -I scripts/make_targets.py  (from PROJECT4/)
"""
import csv
import io
import json
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DL = ROOT / "data/targets/casp16_downloads"
TGT = ROOT / "data/targets"


def read_fasta(path):
    """Return list of (header, sequence) in file order."""
    out, head, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if head is not None:
                out.append((head, "".join(seq)))
            head, seq = line[1:], []
        elif line:
            seq.append(line)
    if head is not None:
        out.append((head, "".join(seq)))
    return out


def entries(path, target_id):
    return [s for h, s in read_fasta(path) if h.split()[0] == target_id]


def smiles_from_tar(tar_path, member, copies=1):
    """Return (name, SMILES) from an official CASP SMILES .tsv. `copies` = expected number of
    rows; multi-copy files (e.g. one ligand per protomer of a homodimer) must list the SAME
    ligand name and SMILES in every row, otherwise this fails rather than guessing."""
    with tarfile.open(tar_path) as tf:
        text = tf.extractfile(member).read().decode()
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    assert len(rows) == copies, f"{member}: expected {copies} ligand row(s), got {len(rows)}"
    assert len({(r["Name"], r["SMILES"]) for r in rows}) == 1, f"{member}: rows are not identical copies"
    return rows[0]["Name"], rows[0]["SMILES"]


def mature_from_IIGG(seq):
    """Chymotrypsin-family serine proteases are activated by cleavage before Ile16
    (chymotrypsin numbering), so the mature enzyme starts at IIGG. The signal
    peptide + activation dipeptide preceding it are absent from the crystal.
    Justified from biology, not from the native file."""
    assert seq.count("IIGG") == 1
    return seq[seq.index("IIGG"):]


T1 = DL / "casp16.T1.seq.txt"
H1 = DL / "casp16.H1.seq.txt"
RDM = DL / "casp16.RDM1.seq.txt"

(t1206,) = entries(T1, "T1206")
h_nb, h_mng = entries(H1, "H0215")  # H1215 is filed as H0215 (same sequences on the H1215 page)
m1211_prot, m1211_rna = entries(RDM, "M1211")
(r1211,) = entries(RDM, "R1211")
m1276_prot, m1276_dna = entries(RDM, "M1276")
(chymase_pre,) = entries(DL / "L1000.targetpage.seq.txt", "L1000")
(catg_pre,) = entries(DL / "L2000.targetpage.seq.txt", "L2000")
(mpro,) = entries(DL / "L4000.targetpage.seq.txt", "L4000")  # mature Mpro, 306 aa; no trimming
assert m1211_rna == r1211

lig1002 = smiles_from_tar(TGT / "L1000.SMILES.tar.gz", "L1000/L1002.tsv")
lig2001 = smiles_from_tar(DL / "L2000.SMILES.tar.gz", "L2000/L2001.tsv")
# L4001 (SARS-CoV-2 Mpro, set L4000): CASP README states the structure is a homodimer; the
# official SMILES file lists 2 identical LIG rows (one per protomer) -> A2 + 2 ligand copies.
lig4001 = smiles_from_tar(TGT / "L4000.SMILES.tar.gz", "L4000/L4001.tsv", copies=2)

# chains: list of {type, ids, sequence|smiles}. ids are the chain ids used in every method's input.
targets = {
    "T1206": {"type": "protein-protein homodimer", "stoich": "A2",
              "chains": [{"type": "protein", "ids": ["A", "B"], "sequence": t1206}]},
    "H1215": {"type": "protein-protein heterodimer", "stoich": "A1B1",
              "chains": [{"type": "protein", "ids": ["A"], "sequence": h_nb, "note": "nanobody"},
                         {"type": "protein", "ids": ["B"], "sequence": h_mng, "note": "mNeonGreen"}]},
    "M1211": {"type": "protein-RNA", "stoich": "A2R1",
              "chains": [{"type": "protein", "ids": ["A", "B"], "sequence": m1211_prot},
                         {"type": "rna", "ids": ["C"], "sequence": m1211_rna}]},
    "R1211": {"type": "RNA", "stoich": "A1",
              "chains": [{"type": "rna", "ids": ["A"], "sequence": r1211}]},
    "M1276": {"type": "protein-DNA", "stoich": "A1B1",
              "chains": [{"type": "protein", "ids": ["A"], "sequence": m1276_prot},
                         {"type": "dna", "ids": ["B"], "sequence": m1276_dna}]},
    "L1002": {"type": "protein-ligand", "stoich": "A1",
              "chains": [{"type": "protein", "ids": ["A"], "sequence": mature_from_IIGG(chymase_pre),
                          "note": "chymase, mature chain from IIGG (precursor 247 aa from CASP L1000 page)"},
                         {"type": "ligand", "ids": ["L"], "smiles": lig1002[1], "ligand_name": lig1002[0]}]},
    "L2001": {"type": "protein-ligand", "stoich": "A1",
              "chains": [{"type": "protein", "ids": ["A"], "sequence": mature_from_IIGG(catg_pre),
                          "note": "cathepsin G, mature chain from IIGG (precursor 255 aa from CASP L2000 page)"},
                         {"type": "ligand", "ids": ["L"], "smiles": lig2001[1], "ligand_name": lig2001[0]}]},
    "L4001": {"type": "protein-ligand", "stoich": "A2 + 2 LIG",
              "chains": [{"type": "protein", "ids": ["A", "B"], "sequence": mpro,
                          "note": "SARS-CoV-2 Mpro homodimer (CASP L4000 README); sequence from CASP L4000 page, 306 aa"},
                         {"type": "ligand", "ids": ["L", "M"], "smiles": lig4001[1], "ligand_name": lig4001[0]}]},
}

out = TGT / "targets.json"
json.dump(targets, open(out, "w"), indent=2)
for k, v in targets.items():
    desc = ", ".join(f"{c['type']}x{len(c['ids'])}"
                     + (f"({len(c['sequence'])})" if "sequence" in c else f"({c['ligand_name']})")
                     for c in v["chains"])
    print(f"{k}: {v['stoich']:5s} {desc}")
print("wrote", out)
