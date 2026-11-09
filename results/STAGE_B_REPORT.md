# Stage B report: seed-1 smoke test (2026-11-09)

All 16 runs (7 AF3, 7 Boltz-2, 2 AF2-Multimer) exited 0. That gives 80 predictions = 16 runs × 5 samples.
Every input type ran: protein–protein, protein–RNA, RNA alone, protein–DNA and protein–ligand (SMILES).
Templates were off everywhere. Ligand input was SMILES only.

Tables (in `results/tables/`):
- `manifest.csv`: one row per prediction
- `timing.csv`: MSA and inference in separate columns, plus `msa_source`
- `confidence.csv`: long format
- `msa_depth.csv`

## Per-run results (seed 1, 5 samples each)

| Target | Method | rc | MSA s (source) | Inference s | Run wall s | Output chains | Confidence (best / range) |
|---|---|---|---|---|---|---|---|
| T1206 | AF3 | 0 | 715 (af3_local_jackhmmer) | 56.1 | 91.9 | A prot 237, B prot 237 | ipTM 0.93, pTM 0.93, pLDDT 94.0 |
| T1206 | Boltz-2 | 0 | 0.8 (server cached) | 85 predict | 321.7* | A 237, B 237 | conf 0.945, ipTM 0.93–0.94, pLDDT 94.9 |
| T1206 | AF2-M | 0 | 170.8 (server cached; ~105 s local startup) | – | 231.7 | A 237, B 237 | ipTM 0.89–0.91, pTM 0.93, pLDDT 93.1 |
| H1215 | AF3 | 0 | 1606 (af3_local) | 51.6 | 70.9 | A 131, B 238 | **ipTM 0.11–0.14**, pTM 0.61, pLDDT 80.8 |
| H1215 | Boltz-2 | 0 | 1.9 (cached) | 79 predict | 316.2* | A 131, B 238 | conf 0.858, **ipTM 0.14–0.86** (samples disagree) |
| H1215 | AF2-M | 0 | 37.4 (MIXED: paired fresh) | – | 429.2 | A 131, B 238 | **ipTM 0.16–0.23**, all models used the maximum 20 recycles |
| M1211 | AF3 | 0 | 1473 (af3_local, protein + RNA) | 52.1 | 73.6 | A 192, B 192, C rna 90 | ipTM 0.66–0.72, pTM 0.79, pLDDT 75.2 |
| M1211 | Boltz-2 | 0 | 0.9 (cached, protein only) | 81 predict | 315.4* | A 192, B 192, C rna 90 | conf 0.69, ipTM 0.33–0.56, pLDDT 75.6 |
| R1211 | AF3 | 0 | 306 (af3_local nhmmer) | 40.7 | 72.0 | A rna 90 | pTM 0.40, pLDDT 67.8 |
| R1211 | Boltz-2 | 0 | none (no RNA MSA in Boltz) | 44 predict | 274.4* | A rna 90 | conf 0.66, pTM 0.45, pLDDT 71.5 |
| M1276 | AF3 | 0 | 657 (af3_local) | 40.6 | 74.5 | A prot 196, B dna 33 | **ipTM 0.18–0.22**, pTM 0.82 |
| M1276 | Boltz-2 | 0 | 1.0 (cached) | 57 predict | 282.4* | A 196, B dna 33 | ipTM 0.84–0.90, conf 0.83 |
| L1002 | AF3 | 0 | 872 (af3_local) | 53.4 | 73.5 | A prot 226, L ligand | ipTM 0.82–0.84, pLDDT 95.2 |
| L1002 | Boltz-2 | 0 | 1.6 (cached) | 67 predict | 289.6* | A 226, L ligand | ligand ipTM 0.96, conf 0.948 |
| L2001 | AF3 | 0 | 870 (af3_local) | 52.0 | 70.2 | A prot 235, L ligand | ipTM 0.91–0.92, pLDDT 92.5 |
| L2001 | Boltz-2 | 0 | 1.7 (cached) | 19 predict | 52.5 | A 235, L ligand | ligand ipTM 0.985, conf 0.951 |

\* **Boltz cold-start:** the first six Boltz runs started together on a freshly built env and paid a one-time Triton/cuEquivariance
kernel compile. The caches `~/.triton/cache` and `~/.cache/cuequivariance-triton` were created at 17:00:04–17:00:31, during those runs.
L2001 ran afterwards with the caches already filled, and its wall time was 52.5 s.
Those six wall times overstate Boltz steady-state cost; Stage C runs will have warm caches.
"predict" is Boltz's own "Predicting DataLoader" time.

Confidence is a model's self-assessment, not accuracy; nothing is scored yet.
Two notable disagreements: H1215 (AF3 ipTM 0.11–0.14 and AF2-M ~0.2, vs Boltz up to 0.86 but with samples disagreeing) and M1276 (AF3 ipTM ~0.2 vs Boltz ~0.87).

## MSA caching
- **AF3:** one CPU MSA pass per target with `--norun_inference`. All 7 seed-1 runs used `<name>_data.json` with `--norun_data_pipeline`.
  **Verified from the logs:** 0 jackhmmer/nhmmer launches and 0 "Running data pipeline" lines in all 7.
  Per-seed cost drops from 5–27 min to ~72 s.
- **Boltz-2:** caching is supported natively. `boltz.main.compute_msa` writes CSVs and the YAML `msa:` field reuses them.
  Runs use no `--use_msa_server`, so a missing MSA is an error rather than a silent re-fetch.
- **ColabFold:** `--msa-only` writes the `.a3m`, which is the input for every seed. The logs show no server calls during prediction.

## MSA depth (`results/tables/msa_depth.csv`)

| Target | AF3 unpaired (paired) | Boltz unpaired (paired) | ColabFold a3m |
|---|---|---|---|
| T1206 | **45 (50)** | **76** | **77** |
| H1215 nanobody / mNeonGreen | 10678 (24703) / **410 (550)** | 11222 / **264** (paired 19) | 12783 |
| M1211 protein / RNA | **548** (9139) / RNA **99** | 948 / no RNA MSA | – |
| R1211 RNA | **99** | none | – |
| M1276 protein / DNA | **513** (325) / – | 1077 | – |
| L1002 | 19065 (50000) | 13293 | – |
| L2001 | 19067 (50000) | 14757 | – |

- **Shallow:** T1206 (<80 sequences in every pipeline; a viral spike), the R1211/M1211 RNA (99), and the H1215 mNeonGreen chain (264–410).
- **Moderate:** M1211 and M1276 protein (500–1000).
- **Paired MSA for H1215 is nearly empty in the MMseqs2 pipelines** (19 paired rows). That is expected for a nanobody–antigen pair with no shared organism.

## MSA server cache status (`msas/msa_server_cache_status.tsv`)
- Boltz: **all 6 protein targets were cache hits** (SUBMIT→COMPLETE with no PENDING/RUNNING).
- ColabFold: T1206 was a cache hit. **H1215 was MIXED**: unpaired cached, paired a fresh search of about 28 s on the server.
- So the ColabFold MSA times are not comparable even with each other, and none of them is comparable to AF3's local search.

## Problems found and how they were handled
1. **Boltz attempt 1 failed** (7/7) with `ModuleNotFoundError: cuequivariance_torch`.
   The env was rebuilt on a CUDA-12 torch with `boltz[cuda]` (decision recorded in CLAUDE.md), and the reruns succeeded.
   The failed outputs are in `predictions/_failed/`.
2. **Maintenance window:** the 24 h CPU time limit blocked scheduling, so it was cut to 10 h. GPU job limits are 2–3 h. Nothing will run into Nov 10 08:00.
3. **AF3 MSA time varies between runs:** T1206 took 446 s in the positive control and 715 s in this pass (same input).
   Shared-filesystem load is a likely cause. Expect ±50% on AF3 MSA times.

## Open item
- L2001 cathepsin G is trimmed only at the N-terminus (235 aa). The mature enzyme may also lack a C-terminal extension.
  This is not acted on; it is flagged for a decision before scoring.

## Stage C cost estimate (not started)
Remaining work: seeds 2–5, giving 28 AF3 runs, 28 Boltz runs, and 2 AF2-M jobs of 4 seeds each.

| Method | Basis | GPU-hours |
|---|---|---|
| AF3 | 28 × ~75 s run wall | **0.6** |
| Boltz-2 | 28 × (predict 19–85 s + ~35 s overhead) with warm caches ≈ 28 × ~95 s; worst case 28 × 320 s cold | **0.75** (worst case 2.5) |
| AF2-M | T1206 4 seeds ≈ 4 × 180 s; H1215 4 × 380 s (20 recycles each) | **0.6** |
| **Total** | | **~2 GPU-h (worst case ~3.7)** |

With 6 L40S cards on `dkhf3-lab-gpu`, wall-clock is about **30–45 min**. No new MSA searches are needed.
Stage D scoring runs on CPU (OpenStructure, USalign, PoseBusters), roughly under an hour for about 400 predictions.
