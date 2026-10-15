# Data dictionary

## SO(3) reference trajectories

`results/reference/so3/` contains CSV and JSON results for `alignment` and
`averaging`, each at 1200 and 1400 decimal digits. No generated files are written
in this directory.

CSV columns:

| Column | Meaning |
|---|---|
| `k` | Iterate index. Indices 0 and 1 are the supplied endpoints. |
| `residual_F` | Frobenius norm of the matrix-valued field. |
| `rho_k` | `log(res[k+1]/res[k])/log(res[k]/res[k-1])`; `--` when undefined. |
| `angular_error` | Angular distance to the known root. |
| `distance_to_X1` | Angular distance to the second initial endpoint. |
| `objective` | Geodesic least-squares objective for averaging; `--` for alignment. |
| `orthogonality_F` | Frobenius norm of `X.T*X-I`. |
| `determinant` | Determinant of the iterate. |

The decimal strings retain 50 significant digits. The JSON files contain the
same rows, input parameters, exact evaluation counts, and step diagnostics.
`diagnostic_log_calls_for_objective` counts logarithms used to evaluate the
averaging objective separately; these do not form part of the secant iteration.
`function_calls` includes both initial endpoints and all transverse evaluations.

Step diagnostics include the secant identity defect, orthonormal-frame defect,
step norm, central-difference step size, and the Frobenius norm of `L+I`.

The pair `X_0,X_1` in the experiment corresponds to the mathematical initialization
`p_{-1},p_0`. Consequently, terminal index 13 represents 12 updates, not 13.

For alignment the exact input vectors are listed in
[semilocal bounds](semilocal_bounds.md). For averaging, the initial displacements
from the common rotation are `(0.050,-0.040,0.030)` and
`(-0.040,0.030,0.020)`. The polar factor of the mean of the six symmetrically
perturbed rotations is exactly `R_star`. The solver itself does not receive the
root or a derivative of the field.

## Bessel input matrices

`experiments/bessel/problem_dense_original.npz` is a fixed, non-pickled NumPy archive.
Its SHA-256 is

```text
140361986c803b86e46d6cb7c4190a5748f61e6a8efec966bcc7d21cabfd21b9
```

| Array | Shape | Meaning |
|---|---|---|
| `C` | `(45,45)` | Linear coefficient matrix. |
| `D` | `(45,45)` | Mixing matrix for scalar Bessel arguments. |
| `A` | `(10,10)` | Fixed skew matrix `-(Rstar-Rstar.T)`. |
| `Rstar` | `(10,10)` | Prescribed root in `SO(10)`. |
| `Xstart` | `(10,10)` | Rotation used only to define the initial direction. |

The kernel uses `g(X)=Gamma(X)+A_upper`, which equals
`Gamma(X)-Gamma(Rstar)`. Here `Gamma(Z)=(Z_ij-Z_ji)_{i<j}` and the basis is
`E_ij=e_i*e_j.T-e_j*e_i.T`, ordered lexicographically. In particular,
`Gamma(hat(z))=2*z`; the geometric derivative map includes this factor.

The common start is `Rstar*exp(0.10*hat(u))`, where `u` is the unit coordinate
vector of `log(Rstar.T*Xstart)`. Only the secant requires the additional point
obtained from a coordinate step of `-0.01*ones(45)/sqrt(45)`.

The field is

\[
f(X)=Cg(X)+0.15\Phi(3\mathbf1+0.4Dg(X))-K,\qquad
K=0.15\Phi(3)\mathbf1,
\]

with componentwise scalar function

\[
\Phi(t)=\sum_{\ell=0}^{99}\frac{1+0.25\cos(0.7\ell)}{100}
J_\ell((\ell+1)t/3).
\]

The root is prescribed for a controlled numerical experiment; its availability
is not a general initialization procedure for unknown-root problems.

## Bessel timing records

`results/reference/bessel/blocks.csv` has 51 rows and four columns:
`block`, `newton`, `secant`, and `ratio`. Times are mean milliseconds per
complete solve over five solves in the corresponding paired block.
`ratio=secant/newton` is calculated separately for each row.

`benchmark.json` contains the same rows, parameters, software environment,
solver counts and residuals, accuracy diagnostics, and summary statistics.
`hardware.json` is the machine record captured before the timing experiment
on 2026-10-04. These records are distributed unchanged.

The reported intervals are the first and third quartiles, not confidence
intervals. The median paired ratio is not the ratio of the two marginal median
times. `verify_reference.py` checks these distinctions directly.

## Scalar data

`results/reference/so2/so2_results.csv` records the scalar equation
`2*sin(theta)+1/4=0`, with initial values `0` and `1/7` at 1800 digits.
Columns are `k`, `theta`, `residual_frobenius`, and `rho_k`.
The first index satisfying the Frobenius tolerance `1e-700` is 13.
Index 14 is a separate diagnostic iterate required to compute `rho_13`.
