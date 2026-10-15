# Reproduction and interpretation

## Environment

The numerical reference environment uses Python 3.13.5, mpmath 1.3.0,
NumPy 2.3.5, and SciPy 1.17.0. `threadpoolctl` controls the number of BLAS
threads. The installation file pins all four third-party dependencies used by
this package. A virtual environment is recommended; the commands do not require
credentials or a network connection after dependencies have been installed.

The reference Bessel timings were collected on an AMD EPYC 9V74 reported by a
KVM virtual machine with five visible logical CPUs, Debian GNU/Linux 13,
Linux 6.18.44, and OpenBLAS 0.3.30. The complete record is in
`results/reference/bessel/hardware.json`. Installing the same Python packages
does not reproduce a hardware environment or guarantee matching wall-clock times.

## Verification levels

`python verify_integrity.py` verifies the bytes of all distributed files
against `SHA256SUMS`. It uses the standard library only. SHA-256 detects changes
relative to this distribution; it is not a proof of scientific correctness or a
digital signature.

`python verify_reference.py --records-only` checks the archived CSV/JSON data,
stored precision agreement, exact rational bounds, metadata agreement, and
statistics of the recorded Bessel timing blocks. It does not execute the solvers.

`python run_all.py --quick` also runs the algebraic tests and the Bessel
solvers. It deliberately skips recomputation of the high-precision trajectories.

`python run_all.py` recomputes the alignment and averaging trajectories at
1200 and 1400 digits, the scalar trajectory at 1800 digits, and the Bessel
solvers. It then compares the new high-precision results with the references.
It does not time the methods.

A new timing experiment must be requested with
`python -m experiments.bessel.benchmark` or `python run_all.py --benchmark`.

## What comparisons require

The high-precision comparison checks all 50 stored significant digits of the
SO(3) residuals, objective values, angular errors, distances, and defined order
indicators. It also checks the update and field-evaluation counts. Very small
orthogonality and determinant errors are treated as roundoff diagnostics; they
are not required to agree bitwise across arithmetic environments.

For the Bessel solvers, the checks require the expected iteration and evaluation
counts, satisfaction of the common tolerance, small secant defects, agreement of
the analytic Jacobian with independent centered directional differences, and
proximity to the known root. Final residuals near double-precision roundoff are
not required to match the recorded values bitwise.

For the recorded timing data, CSV and JSON blocks must agree and the quantiles
must be reproduced. A fresh timing run is not checked against a prescribed speed
ratio, and no failure is inferred merely because it is faster or slower on a
different machine.

## Cost accounting

The SO(3) reduced operator uses four transverse evaluations in dimension three.
Each update needs one further evaluation at the new iterate, for a total
`N_F=2+5*k`. Evaluations of the averaging objective are diagnostics and are
counted separately from field evaluations.

Both Bessel methods start from the same rotation. Fixed matrices, the constant
`K`, and the shared starting rotation are setup, excluded from both timers.
The secant's additional starting point, its exponential, and its field
evaluation are included in its timer. Both timers include residual evaluation,
operator assembly, linear-system solution, and update exponentials.

The Bessel field evaluates 100 modes at 45 component arguments in one vectorized
call. Newton uses one additional such call per analytic Jacobian after reusing
values from the residual. In the recorded trajectories, the total numbers of
full Bessel batches are seven for Newton and six for the secant. `fallbacks`
counts scalar slope substitutions; it is zero in those trajectories.
`bessel_J_batches` counts the additional batches used by Newton's cached
derivative calculation, not calls made by a scalar fallback. The tests instrument
the Bessel function directly to verify that distinction.

The timing procedure uses 20 warmup pairs, then 51 paired blocks, alternating
which method runs first. Each block averages five complete solves per method.
All timed solves use one BLAS thread and wall-clock `perf_counter_ns`, not CPU
process time. Repetitions concern one fixed problem, not a sample of independent
nonlinear problems.

## Failures and diagnostics

Singular linear systems, nonfinite steps, exit from the chosen logarithm branch,
or failure to reach tolerance within the iteration budget must not be silently
replaced by another method. These routines raise errors on these events.
The SO(3) routines are configured for their fixed high-precision problems;
mathematical domain bounds and post-run containment checks are documented in
`semilocal_bounds.md`.

The tests exercise signed Householder frames, secant identities, the fixed-basis
central construction, the second-order midpoint error on a smooth test field,
the `SO(3)` principal-logarithm cut locus, cached Bessel assembly, the scalar
fallback path, and reuse of logarithmic increments. Such tests do not replace
uniform analytic estimates.

The GitHub workflow runs quick verification on pushes and pull requests. A
manually requested workflow can also run the full high-precision suite. It does
not publish releases or send data to Zenodo.
