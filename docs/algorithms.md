# Operators and iteration

## Coordinates and increments

A field takes values in a fixed Lie-algebra coordinate space. For two distinct
points `X,Y` in a logarithm domain, define

\[
v=\log(X^{-1}Y),\quad r=\|v\|,\quad q=v/r,\quad
M=X\exp(v/2),\quad h=\tau r.
\]

In the code, exponential maps take coordinate vectors and logarithms return
coordinate vectors. Matrix multiplication and transpose are the default group
product and inverse for the orthogonal groups. Optional callbacks in the generic
helpers allow other group representations.

The chosen logarithm is inverse to the exponential on an open set `W` containing
zero. An accepted increment `s` must lie in `W`; a bound on its norm is useful only
when it guarantees this inclusion. On `SO(n)`, the test `||s||_2 < pi` is a
sufficient principal-branch condition for the skew-matrix coordinates used here,
not a substitute for every possible problem-specific domain restriction.

## Rank-one correction

For a linear auxiliary operator `A`, define

\[
L=A+\frac{(f(Y)-f(X)-Av)v^T}{r^2}.
\]

Then `Lv=f(Y)-f(X)`. The correction vanishes on the orthogonal complement of `v`.
`codes.lie_divdiff.rank_one_correction` implements this formula. The correction
alone does not guarantee invertibility or convergence.

The fixed-basis central helper forms the columns

\[
D_{M,h}e_j=\frac{f(M\exp(he_j))-f(M\exp(-he_j))}{2h}
\]

and applies the rank-one correction. With endpoint values already available,
this construction uses `2*m` new field evaluations. It does not evaluate the
field at `M`.

## Reduced divided difference

Choose an orthonormal frame `Q=[q_1,...,q_m]` with `q_1=q`. In this frame set

\[
Lq_1=\frac{f(Y)-f(X)}r,\qquad
Lq_j=\frac{f(M\exp(hq_j))-f(M\exp(-hq_j))}{2h},\quad j\ge2.
\]

The matrix in fixed coordinates is `D*Q.T`, where the columns of `D` are the
right-hand sides. Exactly `2*(m-1)` new field evaluations are required. The
`SO(3)` reference data use this construction with `tau=0.1`.

The frame is obtained with a signed Householder transformation. Writing
`s=1` when the first component of `q` is nonnegative and `s=-1` otherwise, let

\[
w=q+se_1,\quad H=I-2ww^T/\|w\|^2,\quad
Q=H\operatorname{diag}(-s,1,\ldots,1).
\]

Since `||w||^2=2*(1+abs(q[0]))`, the denominator is bounded away from zero.
The rule is deterministic but not globally continuous as a function of the
increment direction. The midpoint consistency estimate, rather than global
continuity of the frame, is used for its diagonal limit.

## Cached algorithmic secant

For the Bessel field, cache `f(X)`, `w(X)`, and `Phi(w(X))` at each endpoint. With

\[
\sigma_i=\Phi[w_i(X),w_i(Y)],\quad
B_{X,Y}=C+0.06\operatorname{diag}(\sigma)D,
\]

and `G(Z)s=Gamma(Z*hat(s))`, the matrix is

\[
A=B_{X,Y}\mathcal G((X+Y)/2),\qquad
L=A+\left(\frac{f(Y)-f(X)}r-Aq\right)q^T.
\]

The arithmetic mean `(X+Y)/2` is used only as an ambient matrix. It is not an
iterate and the field is not evaluated there.

A two-point implementation follows this procedure:

```text
Evaluate and store the two endpoint records at X_previous and X_current.
Set v = log(inverse(X_previous) * X_current).
While the current residual is above tolerance:
    Form L from v and the two stored records.
    Stop with failure if L is singular.
    Solve L*s = -f(X_current).
    Form Y = X_current * exp(s).
    Stop with failure if s is nonfinite, s is outside W, or Y is outside the field domain.
    Evaluate and store the record at Y once.
    Replace the endpoint records by the records at X_current and Y.
    Set X_previous = X_current, X_current = Y, and v = s.
```

The assignment `v=s` is justified by
`log(inverse(X_current)*Y)=log(exp(s))=s` on `W`. There is no matrix logarithm
inside the cached update loop. If assembly needs no complete field evaluations,
`k` updates use `k+2` field evaluations. This is a function-evaluation count,
not an assertion that matrix assembly, linear solves, exponentials, or scalar
limits have zero cost.

For the reduced `SO(3)` experiments, each update uses four transverse evaluations
and one evaluation at the new iterate. Thus the total is `2+5*k`.

## Scalar slopes and floating-point arithmetic

An exact scalar divided difference is `(Phi(b)-Phi(a))/(b-a)` for distinct
arguments and `Phi'(a)` for coincident arguments. The Bessel implementation
uses `Phi'((a+b)/2)` when

\[
|b-a|\le 10^{-12}(1+\max(|a|,|b|)).
\]

For distinct arguments this is a numerical approximation to the exact slope.
`fallbacks` counts the number of scalar components receiving this substitution.
Each assembly with such components evaluates two vectorized Bessel calls on
those components. The tests check this branch explicitly. No such substitution
is used along the recorded `SO(10)` trajectory.

Newton reuses the Bessel values stored by the field evaluation through

\[
J_\ell'(z)=\frac{\ell}{z}J_\ell(z)-J_{\ell+1}(z),\qquad
\Phi'(t)=\sum_\ell a_\ell\kappa_\ell J_\ell'(\kappa_\ell t).
\]

The arguments of this fixed experiment stay away from zero. The numerical
kernel is intended for this experiment, not as a general-purpose treatment of
the derivative formula at zero.
