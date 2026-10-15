"""Reduced divided differences in fixed Lie-algebra coordinates.

The routines use only function values and user-supplied group exp/log operations.
mpmath matrices represent Lie-algebra coordinates in a fixed orthonormal basis.
No derivative, prescribed root, or special structure of the field is used here.
"""
from __future__ import annotations
from collections.abc import Callable
import mpmath as mp


def norm(x: mp.matrix) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(z)**2 for z in x))


def adapted_frame(v: mp.matrix) -> mp.matrix:
    """Stable, deterministic Householder completion; not globally continuous.

    Columns form an orthonormal basis and column 0 is v/||v||.
    The sign branch can switch when the first component crosses zero.
    The frame is defined for every nonzero increment.
    """
    v=mp.matrix(v); m=v.rows; r=norm(v)
    if r == 0:
        raise ValueError('Distinct secant points are required.')
    q=v/r
    if m == 1:
        return mp.matrix([[q[0]]])
    sign=mp.mpf(1) if q[0]>=0 else mp.mpf(-1)
    w=mp.matrix(q); w[0]+=sign
    Q=mp.eye(m)-2*w*w.T/(w.T*w)[0]
    for i in range(m):
        Q[i,0]=-sign*Q[i,0]
    return Q


def reduced_divdiff(field: Callable, X, Y, fX: mp.matrix, fY: mp.matrix,
                    exp: Callable, log: Callable, tau=None, *, multiply=None, inverse=None):
    """Return (L,v,Q,h). Exactly 2*(m-1) new calls to field.

    exp accepts coordinate vectors; log returns coordinate vectors.
    Optional multiply/inverse callbacks support other Lie-group realizations.
    Defaults are matrix multiplication and transpose (the SO(n) use case).
    """
    tau=mp.mpf('0.1') if tau is None else mp.mpf(tau)
    if tau<=0:
        raise ValueError('tau must be positive.')
    multiply = (lambda A, B: A*B) if multiply is None else multiply
    inverse = (lambda A: A.T) if inverse is None else inverse
    v=mp.matrix(log(multiply(inverse(X),Y))); r=norm(v); Q=adapted_frame(v)
    m=v.rows; h=tau*r; M=multiply(X,exp(v/2))
    D=mp.zeros(m,m)
    y=mp.matrix(fY)-mp.matrix(fX)
    for i in range(m): D[i,0]=y[i]/r
    for j in range(1,m):
        q=Q[:,j]
        z=(mp.matrix(field(multiply(M,exp(h*q))))-mp.matrix(field(multiply(M,exp(-h*q)))))/(2*h)
        for i in range(m): D[i,j]=z[i]
    return D*Q.T,v,Q,h


def rank_one_correction(A: mp.matrix, v: mp.matrix, delta_f: mp.matrix) -> mp.matrix:
    """Impose L*v=delta_f on a square auxiliary matrix in real coordinates."""
    A, v, delta_f = mp.matrix(A), mp.matrix(v), mp.matrix(delta_f)
    m = v.rows
    if v.cols != 1 or A.rows != m or A.cols != m or delta_f.rows != m or delta_f.cols != 1:
        raise ValueError('The auxiliary matrix and coordinate vectors have incompatible shapes.')
    r = norm(v)
    if not mp.isfinite(r) or r == 0:
        raise ValueError('The increment must be finite and nonzero.')
    return A + (delta_f - A*v)*v.T/(r*r)


def central_divdiff(field: Callable, X, Y, fX: mp.matrix, fY: mp.matrix,
                    exp: Callable, log: Callable, tau=None, *, multiply=None, inverse=None):
    """Fixed-basis central auxiliary with rank-one correction.

    Return (L, v, I, h). The coordinate basis is fixed and orthonormal.
    Exactly 2*m new calls to field are used; endpoint values are supplied.
    This helper is separate from the reduced operator used in the SO(3) data.
    """
    tau = mp.mpf('0.1') if tau is None else mp.mpf(tau)
    if not mp.isfinite(tau) or tau <= 0:
        raise ValueError('tau must be finite and positive.')
    multiply = (lambda A, B: A*B) if multiply is None else multiply
    inverse = (lambda A: A.T) if inverse is None else inverse
    v = mp.matrix(log(multiply(inverse(X), Y)))
    r = norm(v)
    if not mp.isfinite(r) or r == 0:
        raise ValueError('Distinct secant points with a finite increment are required.')
    m = v.rows
    h = tau*r
    M = multiply(X, exp(v/2))
    auxiliary = mp.zeros(m, m)
    basis = mp.eye(m)
    for j in range(m):
        direction = basis[:, j]
        values = (mp.matrix(field(multiply(M, exp(h*direction))))
                  - mp.matrix(field(multiply(M, exp(-h*direction)))))/(2*h)
        auxiliary[:, j] = values
    return rank_one_correction(auxiliary, v, mp.matrix(fY)-mp.matrix(fX)), v, basis, h
