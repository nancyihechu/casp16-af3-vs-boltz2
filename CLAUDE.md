# PROJECT4: CASP16 benchmark, AF3 vs AF2-Multimer vs Boltz-2

Root: `/cluster/pixstor/dkhf3-lab/users/nancy/PROJECT4`
Repo: `nancyihechu/casp16-af3-vs-boltz2` (private). GitHub CLI: `~/bin/gh`
Recon: `../PROJECT4_RECON.md`

## Rules (must hold everywhere)
1. **5 seeds per target, the same values for every method: 1, 2, 3, 4, 5.**
   The same seed integer in different codebases does NOT produce comparable randomness.
   Seeds are there for reproducibility only. The fairness control is the **same number of
   predictions per target per method**, so samples-per-seed must match across methods.
   Verify and report the per-method prediction count before any run at scale.
2. **Templates OFF for every method.**
   - AF3: `"templates": []` on every protein chain, with no MSA fields. The MSA search still runs and template search is skipped (verified in `pipeline.process_protein_chain`).
   - Boltz-2: give no `templates:` block in the YAML.
   - AF2-Multimer (ColabFold): run without `--templates`.
3. **Ligands get the SMILES string only.** `ligand_*.pdb` (and `protein_aligned.pdb`) are
   ground truth. Never give them to a predictor, and never derive SMILES or sequences from them.
4. **Record wall-clock time for every run**, with **MSA search time separate from inference time**.

## Other decisions
- Sequences come only from official CASP16 sources (`data/targets/casp16_downloads/`). No UniProt, and never from native files.
- Each method uses its **own** MSA pipeline. MSAs are not shared. Record the MSA source per run.
- Stoichiometry follows CASP, not the chain count in the native file:

| Slot | Target | CASP stoich | Predict | Methods |
|---|---|---|---|---|
| Protein homodimer | T1206 | A2 | 2 × 237 aa | AF3, AF2-Multimer, Boltz-2 |
| Protein heterodimer | H1215 | A1B1 | nanobody 131 + mNeonGreen 238 | AF3, AF2-Multimer, Boltz-2 |
| Protein–RNA | M1211 | A2R1 | 2 × 192 aa + 90 nt RNA | AF3, Boltz-2 |
| RNA only | R1211 | A1 | 90 nt | AF3, Boltz-2 |
| Protein–DNA | M1276 | A1B1 | 196 aa + 33 nt ssDNA | AF3, Boltz-2 |
| Protein–ligand | L1002 (chymase) | A1 | 1 chain + ligand 293 | AF3, Boltz-2 |
| Protein–ligand | L2001 (cathepsin G) | A1 | 1 chain + ligand 761 | AF3, Boltz-2 |
| Protein–ligand | L4001 (SARS-CoV-2 Mpro, added 2026-11-16) | homodimer + 2 LIG | 2 × 306 aa + 2 copies of the ligand | AF3, Boltz-2 |

  L2001: the native has 2 protein chains, which are 2 copies in the crystal ASU. CASP lists A1,
  so predict 1 chain and score against 1 chain. Never drop a chain the target genuinely has.

### Sequence sources
- T1206: `casp16.T1.seq.txt`
- H1215: `casp16.H1.seq.txt`, entries **`>H0215`**. The H1215 target page shows these same sequences, and no `>H1215` entry exists in the file.
- M1211, R1211, M1276: `casp16.RDM1.seq.txt`
- L1002 / L2001 protein: no entry in `sequences/`. Taken from the CASP target pages L1000 / L2000 (`*.targetpage.seq.txt`), which give the full precursors (247 / 255 aa).
  **Decision (2026-11-09): trim to the mature chain, starting at IIGG.**
  - Why: chymase and cathepsin G are chymotrypsin-family preproproteins. Activation cleaves before Ile16 (chymotrypsin numbering), so the mature enzyme starts at IIGG.
  - The crystal holds only the mature enzyme, and the ligand binds the mature enzyme.
  - Predicted precursor residues would have no counterpart at scoring time. The propeptide also occludes the active site in the zymogen, so a model could place it in the pocket and displace the ligand.
  - The trim is justified from biology, not from the native file. It is implemented in `scripts/make_targets.py:mature_from_IIGG`.
  - L1002 chymase: removed `MLLLPLPLLLFLLCSRAEAGE` (21 aa). Mature 226 aa:
    `IIGGTESKPHSRPYMAYLEIVTSNGPSKFCGGFLIRRNFVLTAAHCAGRSITVTLGAHNITEEEDTWQKLEVIKQFRHPKYNTSTLHHDIMLLKLKEKASLTLAVGTLPFPSQKNFVPPGRMCRVAGWGRTGVLKPGSDTLQEVKLRLMDPQACSHFRDFDHNLQLCVGNPRKTKSAFKGDSGGPLLCAGVAQGIVSYGRSDAKPPAVFTRISHYRPWINQILQAN`
  - L2001 cathepsin G: removed `MQPLLLLLAFLLPTGAEAGE` (20 aa). Mature 235 aa:
    `IIGGRESRPHSRPYMAYLQIQSPAGQSRCGGFLVREDFVLTAAHCWGSNINVTLGAHNIQRRENTQQHITARRAIRHPQYNQRTIQNDIMLLQLSRRVRRNRNVNPVALPRAQEGLRPGTLCTVAGWGRVSMRRGTDTLREVQLRVQRDRQCLRIFGSYDPRRQICVGDRRERKAAFKGDSGGPLLCNNVAHGIVSYGKSSGVPPEVFTRVSSFLPWIRTTMRSFKLLDQMETPL`
    **Known limitation (decided 2026-11-09): the C-terminus is NOT trimmed**, so the construct stays at 235 aa.
    Cathepsin G may also lose a C-terminal extension on maturation. It is kept because:
    - the extension is at the opposite end from the active site and cannot displace the ligand;
    - OpenStructure scores common residues only, so the extra residues drop out of scoring;
    - re-trimming would mean redoing the L2001 MSA and every L2001 run.
- **L4001 (Mpro), added 2026-11-16.** It comes from CASP16 pharma set **L4000** (confirmed from `L4000_mpro_README.txt`: Idorsia Mpro, 25 complexes, "The structure is a homodimer").
  - Selection rule, applied before looking at any prediction:
    - Exclude ligands covalently bound to the protein. The README lists the covalent attachments, so 4 are excluded: **L4003, L4013, L4019, L4023**.
    - Then take the median measured affinity. **No affinity data exists for L4000**: none locally, and the CASP `pharma_ligands/` folder has affinity files only for L1000 and L3000. So the fallback applies, the first by ID: **L4001**. That leaves 21 non-covalent candidates.
  - Sequence: CASP target page L4000 (`casp16_downloads/L4000.targetpage.seq.txt`), 306 aa, mature Mpro. It needs no trimming. The L4001 target page has no sequence of its own, as with L1002 and L2001.
  - Stoichiometry: the README says homodimer, and the official `L4000/L4001.tsv` lists **2 identical LIG rows**, one per protomer. So we predict A2 + 2 ligand copies (model chains A, B, L, M).
    The native has chains A (306 aa) and B (304 aa) plus `ligand_LIG_A_401` and `ligand_LIG_B_401`, 19 heavy atoms each, matching the SMILES.
  - SMILES: `L4000.SMILES.tar.gz` → `L4000/L4001.tsv`: `O=C(Cn1ccnn1)Nc1ccc(Br)c([N+](=O)[O-])c1`.
  - MSAs: one AF3 MSA pass (CPU). The Boltz MSA via `boltz_msa.py L4001` was a server **cache hit** (606 sequences).
- **Canonical input definition: `data/targets/targets.json`**, generated by `scripts/make_targets.py`. Every method's input generator reads this file and nothing else.
- SMILES: `L1000.SMILES.tar.gz` → `L1000/L1002.tsv` (ligand 293), and `casp16_downloads/L2000.SMILES.tar.gz` → `L2000/L2001.tsv` (ligand 761).

## Timing: MSA times are NOT comparable across methods
- **MSA time and inference time are always separate columns. Never sum them into one total.**
- `results/tables/timing.csv` carries an `msa_source` column:
  - `af3_local_jackhmmer`: AF3's own jackhmmer/nhmmer search on CPU nodes against the local AF3 DBs. This is a real search (T1206: 446 s in the positive control, 715 s in the MSA pass).
  - `colabfold_server_cached`: the public ColabFold MMseqs2 server returned an already-computed result.
  - `colabfold_server_fresh`: the server ran a new search.
- **Never put an AF3 MSA time next to a ColabFold-server MSA time as if they were comparable.** The <2 s Boltz MSA calls were server cache hits, not search speed.
- Cache status was checked per request from the server state sequence in the logs: SUBMIT→COMPLETE means a cache hit, and PENDING/RUNNING means a fresh search. The record is in `msas/msa_server_cache_status.tsv`.
  - Boltz: all 6 protein targets were cache hits (R1211 has no Boltz MSA).
  - ColabFold: T1206 was a cache hit; H1215 was MIXED (unpaired cached, paired a fresh ~28 s search).
  - So even the ColabFold MSA times are not comparable with each other. The T1206 ColabFold MSA wall time of 170.8 s was ~105 s of local CPU startup, not search.
- **Boltz cold start:** the first Boltz runs on a fresh env pay a one-time Triton/cuEquivariance kernel compile, cached in
  `~/.triton/cache` and `~/.cache/cuequivariance-triton`. The Stage B seed-1 wall times for 6 of 7 targets (275–320 s) include it;
  L2001 (52.5 s) did not. Treat those 6 as cold-start times. Do not clear these caches between Stage C runs.
- AF3 MSA time varies between runs on the same input (T1206: 446 s, then 715 s), most likely from shared-filesystem load.
- AF3 MSA passes run on `general` CPU nodes. The per-seed runs reuse `<name>_data.json` with `--norun_data_pipeline`, so their `data_pipeline_s` = 0 by construction.

## Tools
- **AF3:** `module load alphafold3/alphafold3_v301_deepmind`. This is genuine DeepMind code, git `v3.0.1-13-ga3cf058`, in the SIF
  `/cluster/software/src/SINGULARITY_IMAGES/alphafold3/ALPHAFOLD_3.0.1/alphafold3_301.sif`.
  It provides aliases (`run_script_a3`, `run_script_a3_unified_mem`), so use `shopt -s expand_aliases` in scripts.
  - **Weights:** `/cluster/VAST/cietest/alphafold3_model_parameters/af3.bin`. That directory belongs to another group and may move.
    **`--model_dir` must be passed explicitly** (the default `~/models` does not exist).
  - DBs: `--db_dir=/cluster/VAST/cietest/alphafold3_database`
  - Input dialect is `alphafold3`, version 2, with entity keys `protein` / `rna` / `dna` / `ligand{smiles}`.
  - Do NOT use the `alphafold3/v1_alphafold3` module (it is a third-party reimplementation).
  - Submit with: `sbatch scripts/run_af3.slurm <json> <outdir>`. It writes `<name>_timing.tsv` (wall, data pipeline, inference).
  - Positive control (2026-11-09, job 19698116, g037): T1206 A2, seed 1. Wall time 542.6 s = MSA 446.0 s + inference 54.2 s
    + ~42 s overhead (startup, featurisation, writing). AF3 gives **5 samples per seed** by default (`--num_diffusion_samples=5`).
- **Boltz-2:** our own env `envs/boltz2` (Python 3.11, conda-forge + pip), **rebuilt 2026-11-09**:
  **boltz 2.2.1 (`boltz[cuda]`), torch 2.14.1+cu126 (CUDA 12.6), cuequivariance_torch 0.13.0,
  cuequivariance_ops_torch_cu12 0.13.0, cuequivariance_ops_cu12 0.13.0.** The L40S nodes have driver 610.57.
  - Install recipe: `pip install torch==2.14.1 --index-url https://download.pytorch.org/whl/cu126`, then `pip install "boltz[cuda]==2.2.1"`.
  - Why rebuilt: the first env (plain `pip install boltz`) pulled torch 2.15.0 built for CUDA 13.2 and had no cuEquivariance kernels. All 7 Boltz smoke runs failed with `ModuleNotFoundError: cuequivariance_torch`.
    The `[cuda]` extra pins CUDA-12 kernel packages, so the env was rebuilt on a CUDA-12 torch, which is the documented configuration.
    We did NOT use `--no_kernels`: wall-clock time is a reported metric, and the slow kernel-free path would be our install choice, not a property of the method.
    The old env is kept at `envs/boltz2_cu13_nokernels_RETIRED`, and the failed outputs are in `predictions/_failed/`.
  Weights cache: `BOLTZ_CACHE=/cluster/pixstor/dkhf3-lab/users/nancy/envs/boltz2_cache`.
  Ignore the cluster module `boltz/v2.0.3`: it is mislabelled and ships 2.1.1.
- **AF2-Multimer:** `envs/colabfold` (`colabfold_batch`). The weights are in `../colabfold_weights/params` (all 5 multimer_v3).
- **Scoring:** `envs/p3_tools` (USalign 20241108, posebusters 0.6.5, rdkit 2026.09.1, Biopython, PyMOL).
  Prefer **OpenStructure** for DockQ, lDDT, LDDT-PLI and **BiSyRMSD** (the official CASP16 ligand metric):
  `singularity exec ../containers/ost.sif ost compare-structures | compare-ligand-structures` (OST 2.12.0).
  Install standalone DockQ only if OST cannot cover a case.

## Exact tool versions used (all runs 2026-11-09, L40S, driver 610.57)
| Tool | Version |
|---|---|
| AlphaFold3 | DeepMind container `alphafold3_301.sif`, git `v3.0.1-13-ga3cf058`, weights `af3.bin` |
| Boltz-2 | boltz 2.2.1 `[cuda]`, torch 2.14.1+cu126, cuequivariance(-torch/-ops-cu12) 0.13.0 |
| AF2-Multimer | ColabFold 1.5.5, alphafold-colabfold 2.3.6, jax/jaxlib 0.4.34, `alphafold2_multimer_v3` params |
| MSA (Boltz, ColabFold) | public ColabFold MMseqs2 server `api.colabfold.com`, 2026-11-09 |
| OpenStructure | 2.12.0 (`containers/ost.sif`): compare-structures, compare-ligand-structures |
| USalign | 20241108 (`envs/p3_tools`) |
| PoseBusters | 0.6.5, `dock` mode |
| RDKit (scoring) | 2026.09.1 |
| Biopython | 1.88 |
| ProLIF (3c) | 2.2.3 |
| PDBFixer / OpenMM (3c protonation) | 1.12.0 / 8.3.1 |
| MDAnalysis (ProLIF dependency) | 2.11.0 |

## Stage D scoring conventions
- Pipeline: `scripts/score_prep.py` → `scripts/score_tasks.py <workers>` → `scripts/score_parse.py`. Raw outputs are in `predictions/scoring_work/out/`.
- Ligand SDFs (both model and reference) take connectivity from coordinates and bond orders from the **official SMILES**. Each one is checked by InChI key against the SMILES.
  Radical flags left over from the conversion are cleared; otherwise PoseBusters would report them as a model error.
- PoseBusters runs in `dock` mode, using the model's own protein. Its redock `rmsd_≤_2Å` check cannot be used, because model and reference sit in different frames. RMSD = OST BiSyRMSD only.
- Ligand success = BiSyRMSD < 2 Å AND every PoseBusters check passed.
- L2001 is scored against reference chain A, the chain closest to ligand 761 (2.6 Å). Chain B is 6.1 Å away, a crystal contact.
- Protein–NA headline = mean DockQ over protein–nucleic-acid interfaces only (`dockq_protein_na_mean`). Native NA chain id = `0`.
- **M1211 RNA scored alone:** `rna_only_tm_score` / `rna_only_lddt` (model chain C vs native chain 0, OST `--lddt --tm-score`, the same as R1211).
  The R1211 native and the M1211 native RNA chain are the same coordinates (USalign TM 1.000, RMSD 0.00, 90 nt, 1904 atoms),
  so `results/tables/rna_alone_vs_complex.csv` isolates the effect of context.
- `results/tables/summary.csv` gives per target × method: mean, median, top-ranked (by the method's own ranking_score / confidence_score) and best-of-25.

## Part 1 diagnostics (2026-11-16) — in `results/tables/diagnostics/`, never in benchmark tables
- **1a M1276:** the native has **one** DNA chain (chain `0`), 8 nt `AGTCAGTC`, residues 10–17, contiguous. That makes it single-stranded.
  We predict one 33-nt ssDNA chain, so the composition MATCHES. Recorded in `diagnostics/m1276_chain_composition.csv`.
- **1b Boltz-2 potentials:** in boltz 2.2.1 the only flag is **`--use_potentials`**. The help text reads "Whether to use potentials for steering. Default is False."
  It sets `BoltzSteeringParams.fk_steering` and `physical_guidance_update`, both of which default to False (`boltz/main.py`).
  **The benchmark runs never passed it, so all benchmark Boltz-2 numbers are potentials-OFF.**
  The diagnostic therefore runs both directions for L1002 and L2001, seed 1, using `scripts/run_boltz2_diag.slurm`:
  - `potentials_off_rerun`: no flag, the literal request. It reproduces the benchmark configuration and doubles as a determinism check.
  - `potentials_on`: with `--use_potentials`.
  Outputs are in `predictions/diagnostics/`, and the table is `diagnostics/boltz2_potentials.csv`.

## Part 3 metrics (2026-11-16)
- **3a:** `fnat` for T1206 and H1215 (single interface, from the OST DockQ output).
- **3b:** `ics_protein_na_mean` and `ilddt_protein_na_mean` for M1211 and M1276, restricted to protein–NA interfaces the same way as `dockq_protein_na_mean`. Per-interface rows are `ics[r1-r2]` and `ilddt[r1-r2]`.
  - ICS = OST `per_interface_ics`.
  - OST only gives global `ilddt` (all inter-chain contacts). So each protein–NA interface is re-scored on its own chain pair with the full-complex chain mapping enforced (`-c`), using `scripts/score_na_interfaces.py`.
  - Check: for M1276 (one interface) the result equals the global `ilddt`.
- **3c ProLIF interaction recovery** (`scripts/score_prolif.py`): ProLIF 2.2.3 with its default interactions (Anionic, CationPi, Cationic, HBAcceptor, HBDonor, Hydrophobic, PiCation, PiStacking, VdWContact), `count=False`.
  - **Hydrogens, identical for native and model:**
    - protein: strip every existing H (CASP natives ship with H), then PDBFixer 1.12.0 `addMissingHydrogens(pH 7.0)`. No missing residues or heavy atoms are added.
    - protein bonds: OpenMM residue-template bonds written as explicit CONECT records, read by RDKit with `proximityBonding=False`. Without this, the crystal close contact in the L2001 native (Arg170 NH2–Arg174 NH1, 1.85 Å) became a bond and made the structure unreadable.
    - ligand: the SMILES-bonded SDF, then RDKit `AddHs(addCoords=True)`.
  - Model residues are mapped to native residues by global sequence alignment per chain pair (L1002 check: 226/226 mapped, same residue type).
  - `plif_recovery` = native (residue, interaction type) pairs reproduced / native pairs.
    For L4001, each model copy is compared with the reference copy OST assigned it to; the result pools both copies, and an unassigned copy counts as zero recovered.
  - **PLIF-WM (PoseBench) is NOT reported.** PoseBench builds fingerprints with PoseCheck on a forked ProLIF (`amorehead/ProLIF@patch-2`), using PoseCheck's own protonation.
    It then takes `scipy.stats.wasserstein_distance` of per-target interaction-type histograms and normalises to [0,1] by the min/max EMD over **all method × target entries** (rank-1 poses only).
    So the value is dataset-relative and depends on their fork. `plif_recovery` is not PLIF-WM and must not be labelled as such.
- **3d/3e** (`scripts/make_summary.py` → `summary.csv`):
  - `ranking_loss` = best_of_25 − top_ranked, oriented ≥ 0 (for RMSD: top − best).
  - `seed_mean_min/max/range` = the per-seed mean of 5 samples, then the range across the 5 seeds.
  - Headline metrics: PP `dockq`, protein–NA `dockq_protein_na_mean`, RNA `tm_score`, ligand `success`.
- `scores.csv` only gains rows: all 5601 pre-Part-3 rows are verified unchanged.
- **L4001 ligand aggregation:** per-copy `bisyrmsd[ligandN]`, `lddt_pli[ligandN]`, `pb_valid[ligandN]` and `success[ligandN]`.
  - Per-prediction `bisyrmsd` / `lddt_pli` = mean over assigned copies.
  - `pb_valid` = every copy valid.
  - **`success` = 1 only if EVERY reference copy is placed < 2 Å and PB-valid.** `success_fraction_copies` is also reported.

## Cluster
- Conda: `module load miniconda3/v26.7.1-1_py314; eval "$(conda shell.bash hook)"`. Always use `-c conda-forge --override-channels`.
- Slurm: `--account=dkhf3-lab --partition=dkhf3-lab-gpu --gres=gpu:L40S:1`. This was used for ALL runs on 2026-11-09.
  **By 2026-11-16 the `dkhf3-lab-gpu` partition no longer existed.** The same L40S nodes (g035–g037) were then only in `gpu_requeue`/`requeue`.
  Runs from 2026-11-16 (Part 1b diagnostics, L4001) use `sbatch -p gpu_requeue` with the scripts' `--gres=gpu:L40S:1` (user decision, 2026-11-16).
  `gpu_requeue` is preemptible (PreemptMode=REQUEUE), so a job's `sacct` Restarts count is recorded, and any restarted run is flagged in the timing table.
- Downloads: run them from the interactive compute node (c092 has outbound network), not from the login nodes.

## Layout
`data/targets` (natives + CASP downloads), `data/inputs` (method inputs), `predictions/`, `results/{tables,figures}`, `scripts/`, `logs/`.
`predictions/` and `logs/` are git-ignored.
