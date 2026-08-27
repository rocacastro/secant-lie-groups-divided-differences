# Manuscript reference outputs

These files record deterministic values used to prepare Tables 1, 3, and 4, plus
the hardware-dependent timing run reported in Table 2.

- The `SO(2)` residual and order CSV files are from the unchanged v1.0.0 scripts.
- `so2_algorithms_efficiency_benchmark_output_paper_reference.txt` stores the
  specific timing measurement printed in the manuscript. Fresh timings will vary.
- The `SO(3)` reference CSV and text files were generated with:
  - `mpmath==1.3.0`
  - 1200 decimal digits
  - `tau = 0.1`
  - tolerance `1e-700`

Use `python verify_reference.py` after generating `results/current/`.
