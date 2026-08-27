# Reproducibility code for the secant method on Lie groups

This repository contains the Python code used to reproduce the numerical tables
in the manuscript **Semilocal Convergence of a Secant-Type Method on Lie Groups via Divided Differences**.

This is the revised `v2.0.0` reproducibility package. The two `SO(3)` experiments
now use derivative-free divided differences based on central function evaluations
at the intrinsic midpoint and an exact rank-one secant correction.

## Requirements

- Python 3.10 or later
- `mpmath==1.3.0`

Install the dependency with:

```bash
python -m pip install -r requirements.txt
```

## Scripts and manuscript tables

| Script | Role in the manuscript | Main generated output |
|---|---|---|
| `so2_algorithm1_matrix_high_precision.py` | Algorithm 1 on `SO(2)`; supports Table 1 | Residuals and observed orders |
| `so2_algorithm2_optimized_high_precision.py` | Algorithm 2 on `SO(2)`; supports Table 1 | Residuals and observed orders |
| `so2_algorithms_efficiency_benchmark.py` | Algorithms 1–2 CPU comparison; supports Table 2 | Residual comparison, order estimates, timing data |
| `so3_alignment_bd_rankone_high_precision.py` | Orientation alignment on `SO(3)`; supports Table 3 | Block-diagonal central-difference core plus exact rank-one secant correction |
| `so3_rotation_averaging_r1_high_precision.py` | Geodesic rotation averaging on `SO(3)`; supports Table 4 | Full central-difference core plus exact rank-one secant correction |
| `so3_rankone_divdiff_utils.py` | Shared high-precision `SO(3)` operations and divided-difference utilities | Imported by the two `SO(3)` scripts |

The `SO(2)` computations use 2500 decimal digits. The revised `SO(3)` computations
use 1200 decimal digits, `tau = 0.1`, and stopping tolerance `1e-700`, matching the
manuscript.

## Running all experiments

```bash
python run_all.py
```

Generated files are written to:

```text
results/current/
```

This directory is created automatically and is ignored by Git, except for its
placeholder file.

## Verifying deterministic results

After running all experiments, compare the deterministic numerical columns with
the archived manuscript-reference outputs:

```bash
python verify_reference.py
```

The verifier checks residuals, order estimates, objective values, determinants,
and orthogonality defects. It intentionally does not check CPU times, because
timings depend on hardware, operating system, Python version, and system load.

## Reference outputs

The directory `results/reference/` contains the outputs used when preparing the
manuscript tables. In particular:

- `so2_bench.txt` records the
  timing run reported in Table 2.
- `so3_align.csv` and `so3_align.txt` record the revised alignment
  experiment used in Table 3.
- `so3_avg.csv` and `so3_avg.txt` record the revised rotation
  averaging experiment used in Table 4.

A fresh run should reproduce the deterministic values to the reported precision.
The CPU-time values in Table 2 are retained as a reference measurement and are not
expected to be hardware independent.

## Main changes from v1.0.0

- The component-wise `SO(3)` divided difference is no longer used.
- The alignment experiment now uses a block-diagonal central-difference core and
  an exact rank-one secant correction.
- The rotation-averaging experiment now uses a full central-difference core and
  an exact rank-one secant correction.
- The revised `SO(3)` experiments reproduce an observed order close to the golden
  ratio.
- `run_all.py`, metadata, documentation, and reference outputs were updated.

See `CHANGELOG.md` for additional details.

## Citation

Citation metadata are provided in `CITATION.cff`. After the GitHub release is
archived by Zenodo, cite the DOI assigned to the corresponding release.

## License

The code is released under the MIT License. See `LICENSE`.
