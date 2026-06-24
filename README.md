# MEANDIR — Python port

A Python conversion of **MEANDIR** (*Mixing Elements ANd Dissolved Isotopes in
Rivers*), the Monte Carlo inversion model for dissolved river chemistry by
**Preston Cosslett Kemeny & Mark Albert Torres**.

> Original model (MATLAB): Kemeny & Torres, *Presentation and applications of
> Mixing Elements ANd Dissolved Isotopes in Rivers (MEANDIR)*, American Journal
> of Science (2021); v1.3 "Alaska" release (2023). The original MATLAB source,
> Excel inputs, and user guide are vendored under [`reference/`](reference/) and
> [`docs/`](docs/) for attribution and validation. All scientific credit for the
> model belongs to the original authors.

## Goal & fidelity statement

This is a faithful port that targets **algorithmic and statistical fidelity**:

- **Identical math.** Every equation, normalization, charge-balance closure,
  fractionation calculation, Cl-critical correction, mass balance, excess-SO₄,
  and R/Z/C/W/Y weathering computation reproduces the MATLAB exactly.
- **Statistical equivalence.** MEANDIR is a Monte Carlo engine whose outputs are
  *distributions* (medians, 5/25/75/95th percentiles). Python results are
  validated to match the MATLAB / published Alaska results within Monte Carlo
  sampling noise.

**What is *not* promised:** bit-for-bit identical numbers on a single run. Two
things make that impossible for any MATLAB→Python port of this model class, and
neither affects the scientific conclusions:

1. **Optimizer.** The `optimize` solvers use MATLAB's `fmincon`; the Python port
   uses SciPy (`SLSQP`/`trust-constr`). Same optimum, different last digits.
2. **RNG.** MATLAB and NumPy have different random number generators, so the
   per-simulation random draws differ even with a fixed seed.

The two linear solvers map almost exactly: `mldivide` → `numpy.linalg.lstsq`,
`lsqnonneg` → `scipy.optimize.nnls` (same Lawson–Hanson algorithm).

## Install

```bash
pip install -e .          # runtime deps: numpy, scipy, pandas, openpyxl
pip install -e ".[test]"  # + pytest
pytest                    # validates the I/O + matrix layer vs the Alaska data
```

## Status

| Phase | Component | State |
|------|-----------|-------|
| 1 | Package scaffold + vendored MATLAB/Excel reference | ✅ done |
| 2 | Foundation: units (`conc2equi`, `delta→ratio`), `go`, river/end-member Excel I/O, river matrix | ✅ done, tested |
| 3 | End-member distribution + sampling engine (`makeEMdistributions`, `PullEndMemberRatios`, scenario parameters, fractionation, Cl-critical) | ⏳ in progress |
| 4 | Monte Carlo inversion loop + solvers (`InvertActiveSimulation`, `CostFunction`, `Master`) | ⏳ |
| 5 | Results/calculations + save (mass balance, fractions, reconstructed obs, excess SO₄, other d34S, end-member values, RZCWY) | ⏳ |
| 6 | Validation against the 5 published Alaska scenarios | ⏳ |
| 7 | Engineer Creek (Yukon) dataset wiring | ⏳ |

## Repository layout

```
meandir/            # the Python package
  constants.py      # variable lists, isotope bookkeeping (Go.m, StandardInput.m)
  go.py             # MEANDIR_Go.m
  datamodel.py      # the `river` container
  river_matrix.py   # MEANDIR_GenerateRiverMatrix.m
  io/               # Excel readers (UserEntries sheets + river data)
  engine/           # (phase 3-4) distributions, sampling, inversion loop
  calc/             # (phase 5) post-inversion calculations
reference/matlab/   # original .m files (read-only reference)
reference/data/     # original Excel inputs (used by the validation tests)
docs/               # original user guide PDF
tests/              # pytest suite
```

## Mapping to the MATLAB files

Each Python module names the MATLAB file(s) it ports in its docstring, so the
two codebases can be cross-checked function by function.
