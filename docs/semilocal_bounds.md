# Exact rational bounds for the SO(3) experiments

## Alignment field and domain

Use the angular metric induced by
`<U,V> = tr(U.T*V)/2` on `so(3)` and the usual hat map from three-dimensional
vectors to skew-symmetric matrices. Thus `||hat(v)||=||v||_2` and its Frobenius
norm is `sqrt(2)*||v||_2`.

Let `R_star` have angle `pi/4` around `(1,1,1)/sqrt(3)`, and define

\[
f(X)=\operatorname{vee}\log(X^TR_*),\qquad
\mathcal U=B(R_*,1/5),\qquad \Omega=B(R_*,1/10).
\]

All logarithms are principal. Choose
`W={hat(v): ||v||_2 < pi}`. The exponential is invertible on `W` onto the
principal-logarithm domain.

The ball `Omega` is geodesically convex. After translating its center to the
identity, its positive-scalar quaternion representatives satisfy
`q_0>cos(1/20)`. A convex combination of two representatives has scalar part
larger than `cos(1/20)` and norm at most one; normalization stays in the same
cap and traces the shorter great-circle arc. Its projection to `SO(3)` stays
inside `Omega`.

For `X,Y` in `Omega`, their angular separation is `r<1/5`. Their midpoint
`M` belongs to `Omega`. With `h=r/10`, all transverse evaluations satisfy

\[
d(M\exp(t\widehat q),R_*)\le d(M,R_*)+|t|<1/10+1/50=3/25
\]

for unit `q` and `|t|<=h`. Thus these evaluations stay in `mathcal U`.

## Uniform operator bound

Write `a=f(Z)` and `theta=||a||_2`. The left-trivialized derivative has matrix

\[
J_Z=-I+\tfrac12\widehat a-\gamma(\theta)\widehat a^2,\qquad
\gamma(\theta)=\frac1{\theta^2}-\frac{1+\cos\theta}{2\theta\sin\theta},
\quad \gamma(0)=1/12.
\]

For `0<=theta<=3/25`, set `x=theta/2`. Since

\[
\gamma(\theta)=\frac{\sin x-x\cos x}{4x^2\sin x},\quad
\sin x-x\cos x=\int_0^x t\sin t\,dt\le x^3/3,\quad
\sin x\ge x(1-x^2/6),
\]

we have `0<=gamma(theta)<=1/[12*(1-theta^2/24)]<1/10`. Consequently,

\[
\|J_Z+I\|_2\le \delta(\theta),\qquad \delta(t)=t/2+t^2/10.
\]

Each column of the reduced divided difference is an average of a directional
derivative along its evaluation arc. Its secant-direction column is bounded
using radius `1/10`; the two transverse columns use radius `3/25`. If `S` is
the reduced matrix, orthonormality gives

\[
\|S+I\|_2^2\le\|S+I\|_F^2
\le\delta(1/10)^2+2\delta(3/25)^2
=\frac{6344217}{625000000}<\left(\frac3{25}\right)^2.
\]

The same bound holds on the diagonal when `S=J_X`. The Banach perturbation
bound and the triangle inequality imply, uniformly for all pairs in `Omega`,

\[
\|S_{X,Y}^{-1}\|_2\le25/22,\qquad
\|S_{X,Y}-S_{P,Q}\|_2\le6/25.
\]

This argument does not assume continuity of the Householder frame away from
the diagonal. Derivatives are used in the analysis, not in the reduced solver.

## Initial data and semilocal constants

The initial rotations are `X_0=R_star*exp(hat(b_0))` and
`X_1=R_star*exp(hat(b_1))`, with exact decimal input vectors

```text
b_0 = (-0.01309040431, -0.03172408338,  0.00427275100)
b_1 = (-0.03234918357,  0.00015621005, -0.00459871328)
```

Their norms satisfy `||b_0||_2<35/1000` and `||b_1||_2<33/1000`. As
`f(X_1)=-b_1`, take

\[
\beta=25/22,\quad \eta=3/80,\quad \omega(s,t)=6/25.
\]

The first correction has norm less than `(25/22)*(33/1000)=3/80`.
With the constant omega, the three scalar factors are

\[
a=b=c=\frac{\beta\omega}{1-\beta\omega}=3/8.
\]

The radius equation and inequalities reduce to

\[
R=\left(1+a+\frac{ab}{1-c}\right)\eta=3/50,\qquad
\beta\omega=3/11<1,\quad c=3/8<1,
\]

and `d(X_1,R_star)+R<93/1000<1/10`. Hence the closed ball of radius `R`
about `X_1` is contained in `Omega`, and the closed algebra ball of radius `R`
is contained in `W`.

`codes/verify_semilocal.py` checks these rational inequalities with
`fractions.Fraction`. It does not infer a uniform bound from numerical samples.
The nonzero cross product of `b_0,b_1` checks that the logarithmic error vectors
are not parallel; it is not used alone to claim that `X_0,X_1` do not commute.

## Rotation averaging

For `R_star=exp(hat(0.45,-0.25,0.35))`, set
`R_i=R_star*exp(hat(xi_i))`, with perturbations
`+/-0.060*e_1`, `+/-0.050*e_2`, and `+/-0.040*e_3` and weights `1/6`.
Opposite perturbations yield a zero of

\[
f(X)=\frac16\sum_{i=1}^6\operatorname{vee}\log(X^TR_i).
\]

Applying the derivative bound separately to the six terms gives

\[
\|df_{R_*}+I\|_2\le\frac16\sum_{i=1}^6
\left(\frac{\|\xi_i\|_2}{2}+\frac{\|\xi_i\|_2^2}{10}\right)
=\frac{7577}{300000}<1.
\]

Thus the derivative at `R_star` is invertible. The objective there is

\[
\mathcal J(R_*)=(0.06^2+0.05^2+0.04^2)/3=77/30000.
\]

These facts are checked exactly. They do not by themselves verify all
semilocal hypotheses for the rotation-averaging starting pair.
