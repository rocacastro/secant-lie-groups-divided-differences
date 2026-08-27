# Manuscript reference outputs

These files record the values used to prepare the numerical tables in the revised manuscript. Short filenames are used deliberately so that browser uploads from Windows preserve the names.

- `so2_a1_res.csv`: SO(2) Algorithm 1 residuals.
- `so2_a1_ord.csv`: SO(2) Algorithm 1 order estimates.
- `so2_a2_res.csv`: SO(2) Algorithm 2 residuals.
- `so2_a2_ord.csv`: SO(2) Algorithm 2 order estimates.
- `so2_cmp_res.csv`: SO(2) residual comparison.
- `so2_cmp_ord.csv`: SO(2) order comparison.
- `so2_bench.txt`: timing output reported for the SO(2) benchmark (hardware dependent).
- `so3_align.csv`, `so3_align.txt`: revised SO(3) alignment experiment.
- `so3_avg.csv`, `so3_avg.txt`: revised SO(3) rotation-averaging experiment.

The SO(3) reference files were generated with 1200 decimal digits, `tau = 0.1`, and tolerance `1e-700`. Run `python run_all.py` and then `python verify_reference.py` to compare newly generated deterministic outputs with these archived values. Timing values are intentionally not checked.
