# Repository update report: v1.0.0 to v2.0.0

## Findings from the v1.0.0 archive

The v1.0.0 archive was internally organized and executable, but its two `SO(3)`
scripts implemented the component-wise divided difference that is no longer used
in the revised manuscript:

```text
so3_secant_example2_restart_high_precision.py
so3_rotation_averaging_componentwise_alg3_high_precision.py
```

Its `run_all.py` and README therefore referred to obsolete `SO(3)` experiments.
The repository metadata also used the GitHub owner `rocastro`, while the public
repository is under `rocacastro`.

## v2.0.0 resolution

The revised package:

- preserves the valid `SO(2)` scripts and timing reference;
- removes the two obsolete `SO(3)` scripts;
- adds the block-diagonal rank-one alignment experiment;
- adds the general rank-one rotation-averaging experiment;
- adds shared high-precision `SO(3)` utilities;
- uses 1200 decimal digits for the revised `SO(3)` tables;
- updates `run_all.py`, README, citation metadata, and Zenodo metadata;
- separates generated outputs from archived manuscript-reference outputs;
- adds an automated deterministic-results verifier;
- corrects the GitHub repository URL to the `rocacastro` account.

## Validation performed

The five scripts were run successfully with `mpmath==1.3.0`. The generated
residuals, order estimates, objective values, determinants, and orthogonality
defects passed the automated comparison against the archived reference CSV files.
Timing values are excluded from deterministic verification.
