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

Agreement is within **~1–2 points for nearly every cell**, worst case ~3–4
(`--success 30`, seed 1). Headline quantities — DIC carbonate/Corg/degassing
split, SO₄ from pyrite (H₂SO₄ production), Cl from precipitation, Na/K from
silicate — all match closely.

### Scenario 2 (CO₂ degassing < 2.5× DIC; the main-text scenario)

| ion | end-member        | published | python | Δ |
|-----|-------------------|----------:|-------:|----:|
| DIC | Carbonate         | 72.7 | 72.2 | −0.5 |
| DIC | Corg oxidation    | 38.5 | 40.6 | +2.1 |
| DIC | Degassing         | −10.7 | −12.4 | −1.7 |
| Ca  | Carbonate         | 75.2 | 76.3 | +1.1 |
| Ca  | Evaporite         | 15.7 | 15.9 | +0.2 |
| Ca  | Silicate          | 4.6 | 1.4 | −3.2 |
| Ca  | Precipitation     | 0.5 | 0.5 | 0.0 |
| Mg  | Carbonate         | 70.4 | 69.1 | −1.3 |
| Mg  | Silicate          | 29.4 | 30.7 | +1.3 |
| Na  | Silicate          | 94.2 | 94.4 | +0.2 |
| K   | Silicate          | 90.6 | 91.1 | +0.5 |
| Cl  | Precipitation     | 100.0 | 100.0 | 0.0 |
| SO₄ | H₂SO₄ production  | 79.6 | 78.5 | −1.1 |
| SO₄ | Evaporite         | 19.8 | 20.9 | +1.1 |
| SO₄ | Precipitation     | 0.6 | 0.6 | 0.0 |

Scenario 1 is comparable (worst deviation ~3.8 points, on the Mg
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

Note: the **ClCritical (cyclic-chloride) correction is *not* involved** here.
All five Alaska scenarios set `PrecProcessing = 'EndMember'` (verified in
`MEANDIR_FindScenarioParameters.m`), so the published inversion treats
precipitation as an ordinary end-member and never invokes
`MEANDIR_ClCriticalCorrection`. The port matches that: ClCritical is implemented
(`engine/clcritical.py`) but gated behind `PrecProcessing == 'ClCrit'` and stays
inert for these scenarios.

Cl, Na, K, SO₄, and the DIC carbon split — the quantities the paper's
conclusions rest on — reproduce to within ~1–2 points, so the port is faithful
to the published inversion.
