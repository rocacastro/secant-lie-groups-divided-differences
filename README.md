# Secant methods on Lie groups

Python code and reference data for solving nonlinear equations on finite-dimensional
Lie groups with two-point secant operators. The experiments cover orientation
alignment and geodesic rotation averaging on `SO(3)`, a Bessel field on `SO(10)`,
and scalar validation on `SO(2)`.

The `SO(3)` experiments use a reduced divided difference obtained from transverse
central function evaluations and the exact secant value in the displacement
direction. The `SO(10)` experiment uses an algorithmic divided difference that
reuses intermediate function values at its two endpoints and compares it with
Newton's method using an analytic Jacobian. A fixed-basis central auxiliary and a
general rank-one correction are also available as tested numerical helpers.

## Installation

Use Python 3.13 and run the commands from the repository root. The reference
computations use Python 3.13.5. Exact dependency versions are listed in
`requirements.txt`; no package installation or compilation of this repository is
required.

```bash
python -m venv .venv
# Linux or macOS:
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

For the high-precision `SO(3)` and scalar commands alone, the only third-party
requirement is `mpmath==1.3.0`, also provided in
`requirements-high-precision.txt`.

## Quick verification

```bash
python run_all.py --quick
```

This checks file integrity, citation metadata, exact rational bounds, reference
CSV/JSON agreement, stored 1200/1400-digit residual agreement, algebraic identities,
and the two Bessel solvers. It recalculates the statistics of the recorded timing
blocks but does **not** time the methods again. It does **not** recompute the
high-precision `SO(3)` or scalar trajectories.

## Full numerical reproduction

```bash
python run_all.py
```

This additionally recomputes both `SO(3)` trajectories at 1200 and 1400 decimal
digits and the scalar trajectory at 1800 digits, then compares the new outputs
with the reference data. High-precision results are stored as decimal strings,
not binary floating-point numbers. This command does not perform a timing
experiment.

Generated files go to `results/generated/` and never replace
`results/reference/`. An alternative output directory can be selected with
`--output-dir`.

| Problem | Command | Reference data |
|---|---|---|
| `SO(3)` alignment | `python -m codes.run_reduced --problem alignment --dps 1200` | `results/reference/so3/alignment_red_1200.*` |
| `SO(3)` rotation averaging | `python -m codes.run_reduced --problem averaging --dps 1200` | `results/reference/so3/averaging_red_1200.*` |
| Exact semilocal bounds | `python codes/verify_semilocal.py` | Constants and derivation in `docs/semilocal_bounds.md` |
| `SO(2)` scalar equation | `python -m supplement.so2_validate` | `results/reference/so2/so2_results.csv` |
| Bessel statistics and solvers | `python -m experiments.bessel.verify` | `results/reference/bessel/` |

For the two `SO(3)` experiments, use `--dps 1400` to reproduce the higher-precision
runs. The reference stopping tolerance is `1e-700`, measured in Frobenius norm,
and the finite-difference parameter is `tau=0.1`.

## Reference results

| Quantity | Alignment | Rotation averaging |
|---|---:|---:|
| Terminal iterate index | 13 | 13 |
| Updates after the two initial evaluations | 12 | 12 |
| Complete field evaluations | 62 | 62 |
| Final Frobenius residual | `1.04916214921e-787` | `6.99054806976e-708` |

In the Bessel experiment the common stopping tolerance is `1e-12` in coordinate
Euclidean norm. Newton uses 3 updates and 4 field evaluations; the algorithmic
secant uses 4 updates and 6 field evaluations. The recorded timing experiment has
51 paired blocks. Its median secant/Newton time ratio is
`0.8807853631014997`. Timing values are specific to the recorded hardware,
software, and instance, and are not universal performance expectations.

A fresh timing experiment is optional:

```bash
python -m experiments.bessel.benchmark
```

It uses 20 warmup pairs, 51 alternating paired blocks, five complete solves per
method in each block, and one BLAS thread in the timed region. It writes separate
files under `results/generated/bessel/`. The input matrices in
`experiments/bessel/problem_dense_original.npz` are fixed; the programs do not generate new
matrices. See [reproduction details](docs/reproduction.md) for the cost model and
accuracy criteria.

## Layout and documentation

- `codes/`: high-precision group operations, divided differences, experiments,
  and exact rational checks.
- `experiments/bessel/`: the fixed Bessel instance, its numerical kernel,
  paired timing procedure, and independent checks.
- `supplement/`: scalar `SO(2)` validation.
- `results/reference/`: numerical trajectories and the recorded timing data.
- `tests/`: algebraic, branch-domain, evaluation-count, and consistency checks.
- `docs/`: [algorithms](docs/algorithms.md),
  [semilocal bounds](docs/semilocal_bounds.md),
  [data dictionary](docs/data_dictionary.md), and
  [reproduction details](docs/reproduction.md).

All automated tests can be run separately:

```bash
python -m unittest discover -s tests -v
python verify_integrity.py
python verify_reference.py
```

After a full run, compare generated high-precision trajectories explicitly:

```bash
python verify_reference.py --generated-dir results/generated
```

The numerical tests supplement the analytic estimates. In particular, sampling
and observed order indicators are not proofs of semilocal convergence or exact
asymptotic order.

## Citation and license

Authorship and software citation information are provided in `CITATION.cff`.
The authors are María del Pilar Astudillo, Rodrigo Castro, and Willy Sierra.
`VERSION`, `CITATION.cff`, and `.zenodo.json` are checked for consistent metadata.

The software is distributed under the MIT License; see `LICENSE`.
