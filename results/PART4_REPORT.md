# Project 4: diagnostics, third ligand target, extended metrics (2026-11-16)

Numbers only. Definitions are in `CLAUDE.md`. Data are in `results/tables/` (`summary.csv`, `scores.csv`, `timing.csv`) and `results/tables/diagnostics/`.

## Run conditions
- **GPU estimate given before starting: ~0.3–0.4 GPU-h.** Used: 14 jobs, **0.37 GPU-h** of allocation (`sacct` elapsed sum: 1349 s).
- The `dkhf3-lab-gpu` partition no longer exists. The new runs (1b, L4001) ran on `gpu_requeue` (user decision), `--gres=gpu:L40S:1`, on nodes g035, g038, g039 and g040.
- Preempted or restarted jobs: **0 of 14** (`sacct -D`). The `timing.csv` columns `partition`, `node` and `slurm_restarts` now cover all 90 runs.
- `scores.csv` only gained rows. All 5601 pre-task rows were verified present and unchanged.

## Part 1: diagnostics
**1a. M1276 chain composition** (`diagnostics/m1276_chain_composition.csv`)

| | Chains | DNA |
|---|---|---|
| Native | A protein (194 res), 0 DNA | **1 DNA chain**, 8 nt `AGTCAGTC`, residues 10–17, single-stranded |
| Predicted | A protein (196 res), B DNA | 1 DNA chain, 33 nt, single-stranded |

The composition matches: the native DNA is positions 10–17 of the predicted 33-mer.

**1b. Boltz-2 physical steering potentials.**
- The flag in the installed boltz 2.2.1 CLI help is `--use_potentials`: "Whether to use potentials for steering. Default is False."
- It sets `fk_steering` and `physical_guidance_update`, both of which default to False.
- **The benchmark runs did not pass it, so all benchmark Boltz-2 results are potentials-off.** No potentials-on benchmark numbers exist.
- So two diagnostic runs were made per target (seed 1, 5 samples): `potentials_off_rerun` (no flag, the requested run) and `potentials_on` (`--use_potentials`). They are scored exactly like the benchmark and are kept only in `diagnostics/boltz2_potentials.csv`.
- `potentials_off_rerun` reproduces the benchmark seed-1 BiSyRMSD and confidence scores to 4–5 decimals.

| Target | Run | BiSyRMSD per sample (0–4) | Mean BiSyRMSD | PoseBusters-valid | Failed PB checks | Success | Run wall s |
|---|---|---|---|---|---|---|---|
| L1002 | benchmark_potentials_off | 0.783, 0.876, 0.988, 0.710, 1.600 | 0.991 | 5/5 | none | 5/5 | 289.6 |
| L1002 | potentials_off_rerun | 0.783, 0.876, 0.988, 0.710, 1.600 | 0.991 | 5/5 | none | 5/5 | 85.5 |
| L1002 | potentials_on | 0.722, 0.954, 1.756, 1.046, 1.902 | 1.276 | 5/5 | none | 5/5 | 118.6 |
| L2001 | benchmark_potentials_off | 1.279, 1.294, 1.320, 1.330, 1.313 | 1.307 | 5/5 | none | 5/5 | 52.5 |
| L2001 | potentials_off_rerun | 1.279, 1.294, 1.320, 1.330, 1.313 | 1.307 | 5/5 | none | 5/5 | 65.5 |
| L2001 | potentials_on | 1.315, 1.296, 1.314, 1.360, 1.337 | 1.324 | 5/5 | none | 5/5 | 88.8 |

## Part 2: third ligand target
- **Set:** CASP16 pharma set **L4000**, SARS-CoV-2 Mpro (Idorsia), 25 complexes. The README states "The structure is a homodimer".
- **Excluded as covalent:** 4 (L4003, L4013, L4019, L4023, per the README's covalent-attachment list).
- **Affinity data:** none exists for L4000 (CASP provides affinity files only for L1000 and L3000), so the rule's fallback applies, the first by ID.
- **Ligand the rule picked: L4001**, `O=C(Cn1ccnn1)Nc1ccc(Br)c([N+](=O)[O-])c1`, 19 heavy atoms.
- **Prediction:** 2 × Mpro (306 aa, CASP L4000 page) + 2 ligand copies. The official SMILES file lists 2 identical rows.
  The native has 2 protein chains (306 and 304 res) and 2 ligand copies, and is scored against both chains and both copies.
- **MSAs:** one AF3 MSA pass (CPU, 691.6 s). The Boltz MSA was a ColabFold-server cache hit.
- **Predictions:** AF3 25, Boltz-2 25. All scored; 0 unassigned ligand copies. Per-prediction BiSyRMSD range 1.874–2.149 Å.
- **Timing** (run on 2026-11-16, labelled `D_added_2026-11-16` in `timing.csv`):

| | Inference per seed | Run wall per seed |
|---|---|---|
| AF3 | 83.8–86.1 s | 104.2–109.0 s |
| Boltz-2 | 81.2–92.8 s | 81.2–92.8 s |

- `success` for L4001 = both copies < 2 Å and PoseBusters-valid. `success_fraction_copies` is reported alongside.

## Part 3: metrics added
- **fnat:** for T1206 and H1215.
- **ICS and interface lDDT:** for M1211 and M1276, restricted to protein–NA interfaces like `dockq_protein_na_mean`.
- **ProLIF interaction recovery:** for L1002, L2001 and L4001. ProLIF 2.2.3, default interactions; H via PDBFixer pH 7.0 for protein and RDKit `AddHs` for ligand, the same for native and model.
- **Ranking loss and per-seed spread:** for all targets.
- **PLIF-WM is not reported.** PoseBench computes it with PoseCheck on a forked ProLIF and normalises the EMD against the min/max over all method × target entries in their evaluation. That cannot be reproduced as defined here. The recovery fraction is not PLIF-WM.

## Part 4: results by target type
Columns:
- *Mean of 25*: over all 25 predictions.
- *Top-ranked*: the method's own top prediction (AF3 / AF2-M `ranking_score`, Boltz-2 `confidence_score`).
- *Best of 25*: the best value among the 25.
- *Ranking loss*: best − top-ranked, oriented ≥ 0.
- *Per-seed mean range*: the mean within each seed, then the min–max across the 5 seeds.
- Lower is better for `irmsd` and `bisyrmsd`.

> **Correction (2026-11-17).** AF3's per-sample `ranking_score` is stored to 2 decimals and ties (all 25 T1206 and
> L4001 predictions read 0.93 / 0.96). The first version of this report broke those ties by order, so AF3's
> "top-ranked" pick was wrong for T1206, R1211, L2001 and L4001. AF3's full-precision `ranking_scores.csv` is now used.
> Only AF3 top-ranked and ranking-loss values changed (means, medians, best-of-25 and seed ranges are unchanged);
> notably AF3's L4001 top-ranked success is 0, not 1.

### Protein–protein

| Target | Metric | Method | Mean of 25 | Top-ranked | Best of 25 | Ranking loss | Per-seed mean range (min–max) |
|---|---|---|---|---|---|---|---|
| T1206 | dockq (headline) | AF3 | 0.936 | 0.939 | 0.946 | 0.007 | 0.929–0.939 (range 0.010) |
| T1206 | dockq (headline) | Boltz-2 | 0.863 | 0.889 | 0.891 | 0.002 | 0.854–0.869 (range 0.014) |
| T1206 | dockq (headline) | AF2-Multimer | 0.918 | 0.930 | 0.952 | 0.022 | 0.848–0.939 (range 0.092) |
| T1206 | fnat | AF3 | 0.899 | 0.918 | 0.934 | 0.016 | 0.880–0.911 (range 0.031) |
| T1206 | fnat | Boltz-2 | 0.858 | 0.913 | 0.913 | 0.000 | 0.848–0.868 (range 0.020) |
| T1206 | fnat | AF2-Multimer | 0.914 | 0.907 | 0.945 | 0.038 | 0.880–0.926 (range 0.046) |
| T1206 | irmsd | AF3 | 0.457 | 0.490 | 0.419 | 0.071 | 0.441–0.466 (range 0.025) |
| T1206 | irmsd | Boltz-2 | 0.854 | 0.800 | 0.737 | 0.063 | 0.806–0.900 (range 0.094) |
| T1206 | irmsd | AF2-Multimer | 0.637 | 0.524 | 0.439 | 0.085 | 0.502–1.108 (range 0.606) |
| T1206 | usalign_tm_mm1 | AF3 | 0.996 | 0.996 | 0.997 | 0.001 | 0.996–0.996 (range 0.000) |
| T1206 | usalign_tm_mm1 | Boltz-2 | 0.983 | 0.983 | 0.987 | 0.004 | 0.982–0.984 (range 0.002) |
| T1206 | usalign_tm_mm1 | AF2-Multimer | 0.991 | 0.993 | 0.995 | 0.002 | 0.985–0.994 (range 0.008) |
| H1215 | dockq (headline) | AF3 | 0.051 | 0.119 | 0.298 | 0.179 | 0.027–0.077 (range 0.050) |
| H1215 | dockq (headline) | Boltz-2 | 0.028 | 0.008 | 0.253 | 0.245 | 0.008–0.078 (range 0.070) |
| H1215 | dockq (headline) | AF2-Multimer | 0.022 | 0.019 | 0.031 | 0.012 | 0.020–0.023 (range 0.003) |
| H1215 | fnat | AF3 | 0.054 | 0.184 | 0.265 | 0.081 | 0.028–0.114 (range 0.086) |
| H1215 | fnat | Boltz-2 | 0.026 | 0.000 | 0.367 | 0.367 | 0.000–0.098 (range 0.098) |
| H1215 | fnat | AF2-Multimer | 0.000 | 0.000 | 0.000 | 0.000 | 0.000–0.000 (range 0.000) |
| H1215 | irmsd | AF3 | 12.123 | 9.801 | 3.179 | 6.622 | 9.816–14.486 (range 4.670) |
| H1215 | irmsd | Boltz-2 | 14.866 | 17.512 | 4.010 | 13.502 | 11.630–17.554 (range 5.924) |
| H1215 | irmsd | AF2-Multimer | 8.814 | 9.524 | 7.416 | 2.108 | 8.594–9.053 (range 0.458) |
| H1215 | usalign_tm_mm1 | AF3 | 0.641 | 0.686 | 0.875 | 0.189 | 0.624–0.666 (range 0.043) |
| H1215 | usalign_tm_mm1 | Boltz-2 | 0.623 | 0.610 | 0.752 | 0.142 | 0.610–0.654 (range 0.044) |
| H1215 | usalign_tm_mm1 | AF2-Multimer | 0.614 | 0.614 | 0.616 | 0.003 | 0.614–0.614 (range 0.001) |

### Protein–nucleic acid (interface metrics restricted to protein–NA interfaces)

| Target | Metric | Method | Mean of 25 | Top-ranked | Best of 25 | Ranking loss | Per-seed mean range (min–max) |
|---|---|---|---|---|---|---|---|
| M1211 | dockq_protein_na_mean (headline) | AF3 | 0.018 | 0.015 | 0.021 | 0.005 | 0.017–0.020 (range 0.003) |
| M1211 | dockq_protein_na_mean (headline) | Boltz-2 | 0.021 | 0.024 | 0.046 | 0.022 | 0.018–0.025 (range 0.007) |
| M1211 | ics_protein_na_mean | AF3 | 0.020 | 0.013 | 0.028 | 0.015 | 0.014–0.024 (range 0.011) |
| M1211 | ics_protein_na_mean | Boltz-2 | 0.014 | 0.000 | 0.118 | 0.118 | 0.000–0.027 (range 0.027) |
| M1211 | ilddt_protein_na_mean | AF3 | 0.045 | 0.045 | 0.048 | 0.003 | 0.044–0.046 (range 0.002) |
| M1211 | ilddt_protein_na_mean | Boltz-2 | 0.025 | 0.019 | 0.062 | 0.043 | 0.018–0.034 (range 0.017) |
| M1211 | lddt | AF3 | 0.762 | 0.752 | 0.773 | 0.021 | 0.759–0.765 (range 0.007) |
| M1211 | lddt | Boltz-2 | 0.743 | 0.726 | 0.758 | 0.032 | 0.735–0.751 (range 0.016) |
| M1211 | tm_score | AF3 | 0.499 | 0.497 | 0.504 | 0.007 | 0.498–0.500 (range 0.003) |
| M1211 | tm_score | Boltz-2 | 0.470 | 0.484 | 0.515 | 0.031 | 0.460–0.487 (range 0.027) |
| M1276 | dockq_protein_na_mean (headline) | AF3 | 0.059 | 0.038 | 0.097 | 0.059 | 0.027–0.069 (range 0.042) |
| M1276 | dockq_protein_na_mean (headline) | Boltz-2 | 0.042 | 0.030 | 0.118 | 0.088 | 0.027–0.051 (range 0.024) |
| M1276 | ics_protein_na_mean | AF3 | 0.088 | 0.000 | 0.184 | 0.184 | 0.000–0.134 (range 0.134) |
| M1276 | ics_protein_na_mean | Boltz-2 | 0.011 | 0.000 | 0.222 | 0.222 | 0.000–0.048 (range 0.048) |
| M1276 | ilddt_protein_na_mean | AF3 | 0.087 | 0.010 | 0.159 | 0.149 | 0.013–0.116 (range 0.104) |
| M1276 | ilddt_protein_na_mean | Boltz-2 | 0.026 | 0.024 | 0.144 | 0.120 | 0.008–0.038 (range 0.030) |
| M1276 | lddt | AF3 | 0.837 | 0.829 | 0.847 | 0.018 | 0.830–0.840 (range 0.010) |
| M1276 | lddt | Boltz-2 | 0.810 | 0.814 | 0.822 | 0.008 | 0.806–0.814 (range 0.009) |
| M1276 | tm_score | AF3 | 0.938 | 0.953 | 0.954 | 0.001 | 0.934–0.940 (range 0.006) |
| M1276 | tm_score | Boltz-2 | 0.940 | 0.945 | 0.974 | 0.029 | 0.927–0.951 (range 0.024) |

### RNA (R1211 alone; M1211 RNA chain scored alone)

| Target | Metric | Method | Mean of 25 | Top-ranked | Best of 25 | Ranking loss | Per-seed mean range (min–max) |
|---|---|---|---|---|---|---|---|
| R1211 | tm_score (headline) | AF3 | 0.436 | 0.432 | 0.486 | 0.054 | 0.424–0.448 (range 0.024) |
| R1211 | tm_score (headline) | Boltz-2 | 0.622 | 0.499 | 0.697 | 0.198 | 0.602–0.637 (range 0.035) |
| R1211 | lddt | AF3 | 0.665 | 0.653 | 0.701 | 0.048 | 0.657–0.678 (range 0.021) |
| R1211 | lddt | Boltz-2 | 0.717 | 0.697 | 0.765 | 0.068 | 0.706–0.723 (range 0.017) |
| M1211 | rna_only_tm_score | AF3 | 0.465 | 0.450 | 0.511 | 0.061 | 0.446–0.493 (range 0.047) |
| M1211 | rna_only_tm_score | Boltz-2 | 0.510 | 0.379 | 0.654 | 0.275 | 0.445–0.600 (range 0.155) |
| M1211 | rna_only_lddt | AF3 | 0.684 | 0.664 | 0.716 | 0.052 | 0.670–0.700 (range 0.030) |
| M1211 | rna_only_lddt | Boltz-2 | 0.684 | 0.679 | 0.753 | 0.074 | 0.679–0.695 (range 0.016) |

### Protein–ligand

| Target | Metric | Method | Mean of 25 | Top-ranked | Best of 25 | Ranking loss | Per-seed mean range (min–max) |
|---|---|---|---|---|---|---|---|
| L1002 | success (headline) | AF3 | 0.760 | 1.000 | 1.000 | 0.000 | 0.400–1.000 (range 0.600) |
| L1002 | success (headline) | Boltz-2 | 0.920 | 1.000 | 1.000 | 0.000 | 0.800–1.000 (range 0.200) |
| L1002 | bisyrmsd | AF3 | 1.528 | 1.267 | 0.546 | 0.721 | 0.950–2.285 (range 1.335) |
| L1002 | bisyrmsd | Boltz-2 | 1.272 | 0.914 | 0.685 | 0.229 | 0.991–1.681 (range 0.690) |
| L1002 | lddt_pli | AF3 | 0.835 | 0.831 | 0.959 | 0.129 | 0.752–0.907 (range 0.155) |
| L1002 | lddt_pli | Boltz-2 | 0.886 | 0.931 | 0.946 | 0.015 | 0.840–0.914 (range 0.075) |
| L1002 | pb_valid | AF3 | 1.000 | 1.000 | 1.000 | 0.000 | 1.000–1.000 (range 0.000) |
| L1002 | pb_valid | Boltz-2 | 1.000 | 1.000 | 1.000 | 0.000 | 1.000–1.000 (range 0.000) |
| L1002 | plif_recovery | AF3 | 0.750 | 0.750 | 0.938 | 0.188 | 0.725–0.787 (range 0.062) |
| L1002 | plif_recovery | Boltz-2 | 0.745 | 0.750 | 0.875 | 0.125 | 0.725–0.775 (range 0.050) |
| L2001 | success (headline) | AF3 | 0.680 | 1.000 | 1.000 | 0.000 | 0.400–1.000 (range 0.600) |
| L2001 | success (headline) | Boltz-2 | 1.000 | 1.000 | 1.000 | 0.000 | 1.000–1.000 (range 0.000) |
| L2001 | bisyrmsd | AF3 | 1.777 | 1.045 | 0.771 | 0.274 | 1.050–2.450 (range 1.399) |
| L2001 | bisyrmsd | Boltz-2 | 1.319 | 1.291 | 1.237 | 0.054 | 1.307–1.331 (range 0.024) |
| L2001 | lddt_pli | AF3 | 0.829 | 0.922 | 0.932 | 0.010 | 0.743–0.922 (range 0.179) |
| L2001 | lddt_pli | Boltz-2 | 0.913 | 0.915 | 0.925 | 0.010 | 0.910–0.917 (range 0.008) |
| L2001 | pb_valid | AF3 | 0.960 | 1.000 | 1.000 | 0.000 | 0.800–1.000 (range 0.200) |
| L2001 | pb_valid | Boltz-2 | 1.000 | 1.000 | 1.000 | 0.000 | 1.000–1.000 (range 0.000) |
| L2001 | plif_recovery | AF3 | 0.669 | 0.667 | 0.833 | 0.167 | 0.533–0.767 (range 0.233) |
| L2001 | plif_recovery | Boltz-2 | 0.773 | 0.778 | 0.889 | 0.111 | 0.733–0.800 (range 0.067) |
| L4001 | success (headline) | AF3 | 0.640 | 0.000 | 1.000 | 1.000 | 0.400–0.800 (range 0.400) |
| L4001 | success (headline) | Boltz-2 | 0.800 | 1.000 | 1.000 | 0.000 | 0.400–1.000 (range 0.600) |
| L4001 | bisyrmsd | AF3 | 1.961 | 1.958 | 1.896 | 0.062 | 1.942–1.977 (range 0.035) |
| L4001 | bisyrmsd | Boltz-2 | 1.937 | 1.932 | 1.874 | 0.057 | 1.908–1.968 (range 0.060) |
| L4001 | lddt_pli | AF3 | 0.760 | 0.769 | 0.789 | 0.020 | 0.755–0.766 (range 0.011) |
| L4001 | lddt_pli | Boltz-2 | 0.766 | 0.769 | 0.784 | 0.015 | 0.761–0.775 (range 0.014) |
| L4001 | pb_valid | AF3 | 1.000 | 1.000 | 1.000 | 0.000 | 1.000–1.000 (range 0.000) |
| L4001 | pb_valid | Boltz-2 | 1.000 | 1.000 | 1.000 | 0.000 | 1.000–1.000 (range 0.000) |
| L4001 | plif_recovery | AF3 | 0.739 | 0.731 | 0.808 | 0.077 | 0.715–0.777 (range 0.061) |
| L4001 | plif_recovery | Boltz-2 | 0.622 | 0.692 | 0.731 | 0.038 | 0.577–0.654 (range 0.077) |
| L4001 | success_fraction_copies | AF3 | 0.760 | 0.500 | 1.000 | 0.500 | 0.700–0.800 (range 0.100) |
| L4001 | success_fraction_copies | Boltz-2 | 0.880 | 1.000 | 1.000 | 0.000 | 0.700–1.000 (range 0.300) |
