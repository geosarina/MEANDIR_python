# Validation raw outputs

This folder holds machine-readable (CSV) outputs of the validation tests that
compare the Python port of MEANDIR against the published MATLAB results in
Kemeny et al. (2023, *Global Biogeochemical Cycles* 2022GB007644). Each CSV is
described below: what it tests, the exact run configuration (samples /
simulations / successes / solver settings), and a definition of every column.

The narrative interpretation of these numbers lives in
[`../docs/validation.md`](../docs/validation.md). The CSVs here are the raw
tables behind that narrative, regenerable with the scripts noted per file.

---

## `tableS2_comparison_200successes.csv`

**What it tests.** The headline validation: the inversion-constrained
**fractional contribution of each end-member to each dissolved ion**, compared
cell-by-cell against **Table S2** of the paper's supplement ("The mean of median
end-members fractional contribution to river dissolved load"). This is the
quantity the paper reports and on which its conclusions rest.

**Reference values.** The `published_pct` column is taken from Table S2,
extracted from the supplement PDF by coordinate-based parsing
(`scripts/extract_tableS2.py`) and cross-checked against the values that are
legible in the raw text. Table S2 reports three scenarios × three sample groups.

**How the Python values were produced.**

| setting | value |
| --- | --- |
| script | `scripts/validate_groups.py` (export via `scripts/export_validation_csv.py`) |
| command | `python -m scripts.validate_groups --success 200 --max-iter 15000` |
| scenarios | S1 = no CO₂ degassing; S2 = degassing < 2.5× DIC (main-text); S3 = degassing < 25× DIC |
| sample groups | Mainstem (26 `river` samples), Slough (4 `slough` samples), All (the two together, n = 30) |
| sample selection | functional samples (complete Ca, Mg, Na, K, Cl, SO₄, DIC, δ³⁴S, δ¹³C) classified by the `ExtraField1` river-type label; lake/rhizon/huslia samples excluded, matching the paper's grouping |
| successes per sample | **200** (matching the simulation count used in the paper) |
| max iterations / sample | 15 000 |
| max zero-hits / sample | 4 000 (give up on a sample after this many consecutive non-successful draws) |
| RNG seed | 1 (deterministic / reproducible) |
| solver | `mldivide_optimize` → SciPy SLSQP, `ftol = 1e-6` (matches MATLAB `fmincon`'s `optimset` default) |
| statistic | **mean over samples of the per-sample median** fractional contribution (the paper's "mean of median"), in percent |

**Row count.** 162 = 3 groups × 3 scenarios × 18 (ion, end-member) cells. The
Degassing rows for scenario S1 are blank (`published`/`python`/`delta` empty)
because S1 has no degassing end-member.

**Columns.**

| column | meaning |
| --- | --- |
| `group` | sample group: `All samples`, `Mainstem`, or `Slough` |
| `n_samples` | number of river samples in the group (30 / 26 / 4) |
| `scenario` | `S1`, `S2`, or `S3` |
| `degassing_constraint` | the CO₂-degassing bound for the scenario: `none`, `<2.5x DIC`, `<25x DIC` |
| `ion` | dissolved species the contribution is *to* (DIC, Ca, Mg, Na, K, Cl, SO4) |
| `end_member` | source end-member the contribution is *from* (Carbonate, Corg oxidation, Degassing, Evaporite, Silicate, Precipitation, H2SO4 production) |
| `published_pct` | Table S2 value (%), the MATLAB result; blank where not reported (S1 degassing) |
| `python_pct` | this port's value (%), mean-of-median over the group's samples |
| `delta_pct` | `python_pct − published_pct` (percentage points); blank where no published value |

**How to read it.** Data-constrained quantities (Cl, Na, K from silicate, the
DIC carbonate/Corg/degassing split) have `|delta_pct| ≤ 1` across all groups and
scenarios. The larger deviations are confined to the under-constrained splits
(carbonate↔silicate Ca/Mg, and SO₄ pyrite↔evaporite in the degassing scenarios)
and to the Slough degassing median; `docs/validation.md` explains each as a
small fraction (~5%) of the inversion's own uncertainty band on that parameter.

**Caveats.** Values are a single deterministic realization (seed = 1) at 200
successes; the under-constrained splits carry ~1–2 points of run-to-run Monte
Carlo variance, and the Slough group (n = 4) is the noisiest.
