# Reproducibility code for the secant method on Lie groups

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21092606.svg)](https://doi.org/10.5281/zenodo.21092606)
This directory contains the Python scripts used to reproduce the numerical tables
in the manuscript **Semilocal Convergence of a Secant-Type Method on Lie Groups via Divided Differences**.

## Requirements

The scripts use only the Python standard library and `mpmath`.

```bash
pip install -r requirements.txt
```

The computations are high precision. Running the full benchmark may take a few minutes depending on the machine.

## Scripts and manuscript tables

| Script | Purpose in the manuscript | Main output |
|---|---|---|
| `so2_algorithm1_matrix_high_precision.py` | Algorithm 1 on \(SO(2)\) | Residuals and observed orders for Algorithm 1 |
| `so2_algorithm2_optimized_high_precision.py` | Algorithm 2 on \(SO(2)\) | Residuals and observed orders for Algorithm 2 |
| `so2_algorithms_efficiency_benchmark.py` | Comparison of Algorithms 1 and 2 | Residual comparison, order estimates, CPU-time benchmark |
| `so3_secant_example2_restart_high_precision.py` | Algorithm 3 on \(SO(3)\), orientation alignment | Residual table for the alignment example |
| `so3_rotation_averaging_componentwise_alg3_high_precision.py` | Algorithm 3 on \(SO(3)\), geodesic rotation averaging | Residual, objective and orthogonality-defect table |

The first three scripts support the `SO(2)` comparison tables. The last two scripts support the two `SO(3)` examples.

## Running all scripts

```bash
python run_all.py
```

All generated files are written to:

```text
results/current/
```

This directory is created automatically next to the scripts, regardless of the shell working directory.

## Notes on timing values

The residuals and order estimates are deterministic. CPU and wall-clock timing values depend on hardware, operating system, Python version and system load. For this reason, a fresh run may not reproduce the exact CPU seconds printed in the manuscript, although it should reproduce the same qualitative conclusion: Algorithm 2 has a lower average runtime than Algorithm 1 in the `SO(2)` test.

## Reference output

The directory `results/reference/` contains a text output from the timing run used when preparing the manuscript table. New runs write their own outputs to `results/current/`.

## Citation

If you use this code, please cite the archived Zenodo release:

Rodrigo Castro, María del Pilar Astudillo, and Willy Sierra.  
*Reproducibility code for the secant method on Lie groups via divided differences*.  
Zenodo, 2026. DOI: 10.5281/zenodo.21092606.
