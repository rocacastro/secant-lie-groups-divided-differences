# Changelog

## v2.0.0 — 2026-08-27

### Changed

- Replaced the two obsolete component-wise `SO(3)` implementations used in
  v1.0.0.
- Added `so3_alignment_bd_rankone_high_precision.py`, using a block-diagonal
  central-difference core and an exact rank-one secant correction.
- Added `so3_rotation_averaging_r1_high_precision.py`, using the general
  central-difference core and an exact rank-one secant correction.
- Added `so3_rankone_divdiff_utils.py` with shared high-precision `SO(3)`
  exponential, logarithm, projection, and divided-difference operations.
- Set the revised `SO(3)` computations to 1200 decimal digits to match the
  manuscript tables.
- Updated `run_all.py`, the README, citation metadata, and Zenodo metadata.
- Added deterministic reference outputs and `verify_reference.py`.

### Retained

- The three `SO(2)` scripts and the timing reference used in Tables 1–2 remain
  unchanged.

## v1.0.0 — 2026-06-30

- Initial public reproducibility package.
