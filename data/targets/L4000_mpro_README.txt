The Mpro dataset, generously provided by Idorsia Pharmaceuticals Ltd (Allschwil, CH), contains in total 25 protein-ligand(s) complexes. The structure is a homodimer.

Tasks:
- 25 protein-ligand(s) complexes: ONLY Pose (P) prediction.

The stereochemistry is known and is provided for all the ligands, however the protonation might not be optimal and would need to be determined. Some protein-ligand(s) complexes have more than one ligand present in the active site (20 protein-ligand(s) complexes in total). Molecules located nearby (< 4.5 A) the "main" ligand(s) and interacting with specific or non-specific interactions were also included in the SMILES files (.tsv).

Some ligands are covalently bound to a cysteine residue. Here the attachment points for each covalent inhibitor:
- L4003, LIG, SMARTS: [cX3][nX2][cX3][nX2] (first atom)
- L4013, LIG, SMARTS: [C](=[N])[c] (first atom)
- L4019, LIG, SMARTS: [C](=[N])[c] (first atom)
- L4023, LIG, SMARTS: [C](=[N])[c] (first atom)

In L4007, Cys145 is chemically modified to S-Hydroperoxycysteine (https://pubchem.ncbi.nlm.nih.gov/compound/S-Hydroperoxycysteine).

Reference:
@article {Hazemann2024.02.12.579977,
	author = {Julien Hazemann and Thierry Kimmerlin and Roland Lange and Aengus Mac Sweeney and Geoffroy Bourquin and Daniel Ritz and Paul Czodrowski},
	title = {Identification of SARS-CoV-2 Mpro inhibitors through deep reinforcement learning for de novo drug design and computational chemistry approaches},
	year = {2024},
	doi = {10.1101/2024.02.12.579977},
	URL = {https://www.biorxiv.org/content/early/2024/02/13/2024.02.12.579977},
	eprint = {https://www.biorxiv.org/content/early/2024/02/13/2024.02.12.579977.full.pdf},
	journal = {bioRxiv}
}
