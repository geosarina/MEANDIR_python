# Validation against Kemeny et al. (2023)

The Python port is validated against **Table S2** of the supplement to
Kemeny et al. (2023) (*2022GB007644*), "Inversion-constrained end-member
contributions to river dissolved load." Table S2 reports, per scenario, the
**mean over rivers of the median fractional contribution** of each end-member to
each dissolved ion (the "All samples" column). This maps exactly onto
`mean_over_samples( summary.fraction[ion][em]["median"] )` in this port.

Reproduce with:

```
python -m scripts.validate_against_tableS2 --success 30 --max-iter 12000
```

## Result (published % vs Python %, "All samples")

Agreement is within **~1–2 points for nearly every cell**, worst case ~3–4 on
the carbonate/silicate Ca-Mg split. Headline quantities — DIC carbonate/Corg/
degassing split, SO₄ from pyrite (H₂SO₄ production), Cl from precipitation, Na/K
from silicate — all match closely.

### Scenario 2 (CO₂ degassing < 2.5× DIC; the main-text scenario), 200 successes/sample

| ion | end-member        | published | python | Δ |
|-----|-------------------|----------:|-------:|----:|
| DIC | Carbonate         | 72.7 | 72.0 | −0.7 |
| DIC | Corg oxidation    | 38.5 | 40.0 | +1.5 |
| DIC | Degassing         | −10.7 | −11.5 | −0.8 |
| Ca  | Carbonate         | 75.2 | 76.3 | +1.1 |
| Ca  | Evaporite         | 15.7 | 16.4 | +0.7 |
| Ca  | Silicate          | 4.6 | 1.0 | −3.6 |
| Ca  | Precipitation     | 0.5 | 0.5 | 0.0 |
| Mg  | Carbonate         | 70.4 | 69.1 | −1.3 |
| Mg  | Silicate          | 29.4 | 30.7 | +1.3 |
| Na  | Silicate          | 94.2 | 94.2 | 0.0 |
| K   | Silicate          | 90.6 | 91.1 | +0.5 |
| Cl  | Precipitation     | 100.0 | 100.0 | 0.0 |
| SO₄ | H₂SO₄ production  | 79.6 | 78.0 | −1.6 |
| SO₄ | Evaporite         | 19.8 | 21.3 | +1.5 |
| SO₄ | Precipitation     | 0.6 | 0.6 | 0.0 |

Scenario 1 is comparable (worst deviation ~3.4 points, on the Mg
carbonate/silicate split).

## Interpretation of the residual deviations

The few larger gaps are always a **compensating carbonate ↔ silicate
trade-off** on Ca/Mg (e.g. scenario 1 Mg: carbonate −3.8 / silicate +3.8,
summing to zero). That axis is the most degenerate in the inversion (carbonate
and silicate both source Ca and Mg), so it is the most sensitive to:

1. **Monte Carlo sampling** — this run uses 30 successful simulations per
   sample; the paper uses 200. The "mean of median" carries sampling noise that
   shrinks with more successes.
2. **Optimizer differences** — SciPy SLSQP here vs MATLAB `fmincon`.

### Convergence test at 200 successes/sample (matching the paper)

Re-running at 200 successes (vs 30) resolves whether the gaps are Monte Carlo
noise. They are **not**:

- Well-constrained quantities **tightened** toward the published values — Na
  silicate 94.2 → 94.2 (exact), Cl 100.0 exact, K within 0.5, SO₄ within ~1.5,
  and DIC degassing moved from −12.4 (N=30) to −11.5 (N=200) vs −10.7 published.
- The **carbonate ↔ silicate Ca/Mg split did not converge** — it held at
  ~3–3.6 points (scenario 1 Mg: −3.4/+3.3; scenario 2 Ca silicate 4.6 → 1.0).

So the residual is **systematic, not sampling noise**. It is confined to the one
genuinely under-determined axis of these scenarios: nothing in the observation
set (only δ³⁴S and δ¹³C isotopes — no Ca or Mg isotopes) directly constrains
how Ca and Mg partition between carbonate and silicate, so that split is fixed
mainly by the end-member priors and the optimizer's behaviour in a flat region
of the cost surface. The most likely contributors are the **optimizer**
(SciPy SLSQP vs MATLAB `fmincon`) and the **linear-solve initial condition**
(`numpy.linalg.lstsq`'s minimum-norm solution vs MATLAB `mldivide`'s basic
solution), which can settle into slightly different points along that
degeneracy. This cannot be removed by more Monte Carlo sampling and would
require replicating MATLAB's exact optimizer to eliminate.

Note: the **ClCritical (cyclic-chloride) correction is *not* involved** here.

### Initial-condition experiment (mldivide vs minimum-norm)

We tested whether the residual comes from the linear-solve initial condition fed
to the optimizer. Scenario 2's system is **underdetermined** (9 observations,
10 end-members), so MATLAB's `A\b` returns a *basic* solution (some end-members
exactly zero) while `numpy.linalg.lstsq` returns the *minimum-norm* solution.
`scripts/ab_initial_condition.py` runs both as the optimizer's start and scores
each against Table S2:

- **Scenario 1** (square 9×9): the two are *identical* (for a square full-rank
  system `mldivide` == `lstsq`), so the initial condition cannot explain its
  residual.
- **Scenario 2** (9×10): the literal basic solution is **far worse** — it zeros
  carbonate at the start and the optimizer settles at a wrong minimum, blowing
  the worst Table S2 deviation from **3.4 to 42 points** (DIC carbonate 72→40%).

So MATLAB's *effective* behaviour matches the **minimum-norm** solution, which is
what this port uses. The ~3-point carbonate/silicate Ca-Mg residual is therefore
**not** an initial-condition artifact; it is intrinsic to the optimizer
(SLSQP vs `fmincon`) and end-member sampling on the one under-determined axis,
and cannot be reduced by changing the solve or adding Monte Carlo samples.
All five Alaska scenarios set `PrecProcessing = 'EndMember'` (verified in
`MEANDIR_FindScenarioParameters.m`), so the published inversion treats
precipitation as an ordinary end-member and never invokes
`MEANDIR_ClCriticalCorrection`. The port matches that: ClCritical is implemented
(`engine/clcritical.py`) but gated behind `PrecProcessing == 'ClCrit'` and stays
inert for these scenarios.

Cl, Na, K, SO₄, and the DIC carbon split — the quantities the paper's
conclusions rest on — reproduce to within ~1–2 points, so the port is faithful
to the published inversion.
