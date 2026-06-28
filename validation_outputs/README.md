# Validation raw outputs

This folder holds machine-readable (CSV) outputs of the validation tests that
compare the Python port of MEANDIR against the published MATLAB results in
Kemeny et al. (2023, *Global Biogeochemical Cycles* 2022GB007644). Each test is
exported at up to four granularities — **per-simulation** (raw; both the focused
Table S2 cells and the complete end-member × observation matrix), **per-sample**
(percentile summary), and **summary** (mean-of-median vs published) — and every
file is described below: what it tests, the exact run configuration (samples /
simulations / successes / solver settings), and a definition of every column.

The narrative interpretation lives in
[`../docs/validation.md`](../docs/validation.md); the CSVs here are the raw
tables behind it. Regenerate with the scripts noted per test.

---

## Table S2 fractional-contribution comparison (200 successes/sample)

**What it tests.** The headline validation: the inversion-constrained
**fractional contribution of each end-member to each dissolved ion**, compared
against **Table S2** of the paper's supplement ("The mean of median end-members
fractional contribution to river dissolved load") — the quantity the paper's
conclusions rest on.

**Reference values.** The published values are Table S2, extracted from the
supplement PDF by coordinate-based parsing (`scripts/extract_tableS2.py`) and
cross-checked against the legible raw text.

**Run configuration (all three files below).**

| setting | value |
| --- | --- |
| capture/export script | `scripts/export_tableS2_raw.py` |
| command | `python -m scripts.export_tableS2_raw --success 200 --max-iter 15000` |
| scenarios | S1 = no CO₂ degassing; S2 = degassing < 2.5× DIC (main-text); S3 = degassing < 25× DIC |
| sample groups | Mainstem (26 `river` samples), Slough (4 `slough` samples); "All" = the two together (n = 30) |
| sample selection | functional samples (complete Ca, Mg, Na, K, Cl, SO₄, DIC, δ³⁴S, δ¹³C), classified by the `ExtraField1` river-type label; lake/rhizon/huslia excluded, matching the paper |
| successes per sample | **200** (the paper's simulation count) |
| max iterations / sample | 15 000 |
| max zero-hits / sample | 4 000 |
| RNG seed | 1 (deterministic / reproducible) |
| solver | `mldivide_optimize` → SciPy SLSQP, `ftol = 1e-6` (matches MATLAB `fmincon`'s `optimset` default) |
| cells | the 18 (ion, end-member) pairs Table S2 reports; the **Degassing** cell is absent in S1 (no degassing end-member) and is omitted there |

### `tableS2_persimulation_s200.csv` — raw, one row per simulation × cell

Every successful simulation's fractional contribution, in long format
(≈ 320k rows: 3 scenarios × 30 samples × ≤200 successes × 18 cells, minus S1
degassing and any sample reaching < 200 successes).

| column | meaning |
| --- | --- |
| `scenario` | `S1` / `S2` / `S3` |
| `degassing_constraint` | `none` / `<2.5x DIC` / `<25x DIC` |
| `group` | `Mainstem` or `Slough` (the river-type group; "All" = both combined) |
| `sample_index` | 0-based row index of the sample in the river spreadsheet |
| `river_name` | sample name (e.g. `KY18-AKW-01`) |
| `sim_index` | 0-based index of the successful simulation for that sample (0…n−1) |
| `ion` | dissolved species the contribution is *to* (DIC, Ca, Mg, Na, K, Cl, SO4) |
| `end_member` | source end-member it is *from* (Carbonate, Corg oxidation, Degassing, Evaporite, Silicate, Precipitation, H2SO4 production) |
| `contribution_pct` | the fractional contribution for that single simulation (%) |

### `full_persimulation_s200.csv` — raw, the complete 9×10 matrix per simulation

The same per-simulation data as above but **wide** and **complete**: one row per
successful simulation, with a column for the contribution of every end-member to
every observation (the full inversion result). ≈ 18k rows.

| column | meaning |
| --- | --- |
| `scenario`, `degassing_constraint`, `group`, `sample_index`, `river_name`, `sim_index` | as above |
| `<ion>.<end_member>` | 90 columns, named `ion.end_member`, giving the fractional contribution (%) of that end-member to that observation for the single simulation. Observations: Ca, Mg, Na, K, Cl, SO4, DIC, d34S, d13C. End-members: prec, carb, slct_Ca, slct_Mg, slct_Na, slct_K, pyri, evap, corg, degas. A value of `0.0` means that end-member sources none of that observation; a **blank** means the end-member is absent in that scenario (`*.degas` in S1). |

This is the rawest output: the 18 Table S2 cells above are a subset of these
columns (e.g. `Mg.carb`, `SO4.pyri`). The d34S/d13C columns are the contributions
to the isotope-product rows used internally by the inversion.

### `tableS2_persample_s200.csv` — per-sample percentile summary

The per-simulation values reduced to a distribution per (sample, cell):
≈ 1 600 rows (3 × 30 × 18, minus S1 degassing).

| column | meaning |
| --- | --- |
| `scenario`, `degassing_constraint`, `group`, `sample_index`, `river_name`, `ion`, `end_member` | as above |
| `n_successes` | number of successful simulations for that sample/scenario |
| `median_pct` | median contribution across that sample's simulations (%) |
| `p05_pct`, `p25_pct`, `p75_pct`, `p95_pct` | 5th/25th/75th/95th percentiles (%) — the per-sample uncertainty band |

### `tableS2_summary_s200.csv` — mean-of-median vs published

The headline comparison: 162 rows (3 groups × 3 scenarios × 18 cells). The
`python_pct` is the **mean over the group's samples of each sample's median**
(the paper's "mean of median"); Degassing in S1 is blank.

| column | meaning |
| --- | --- |
| `group` | `All samples`, `Mainstem`, or `Slough` |
| `n_samples` | samples in the group (30 / 26 / 4) |
| `scenario`, `degassing_constraint`, `ion`, `end_member` | as above |
| `published_pct` | Table S2 value (MATLAB result); blank for S1 degassing |
| `python_pct` | this port's mean-of-median (%) |
| `delta_pct` | `python_pct − published_pct` (percentage points) |

**How to read it.** Data-constrained quantities (Cl, Na, K from silicate, the
DIC carbonate/Corg/degassing split) have `|delta_pct| ≤ 1`. The larger
deviations are confined to the under-constrained splits (carbonate↔silicate
Ca/Mg; SO₄ pyrite↔evaporite in the degassing scenarios) and the Slough
degassing median; `docs/validation.md` shows each is a small fraction (~5%) of
the inversion's own uncertainty band — visible directly here in the wide
`p05`→`p95` spread of the per-sample file.

**Caveats.** Single deterministic realization (seed = 1) at 200 successes; the
under-constrained splits carry ~1–2 points of run-to-run Monte Carlo variance,
and the Slough group (n = 4) is the noisiest.
