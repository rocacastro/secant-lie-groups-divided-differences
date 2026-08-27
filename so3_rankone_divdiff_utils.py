"""Shared high-precision SO(3) utilities and rank-one divided differences.

The numerical methods in this module use only evaluations of F, together with
SO(3) exponential/logarithm maps.  No derivative of F is evaluated.
"""
from __future__ import annotations

import mpmath as mp

from _repo_utils import results_dir, write_csv, write_text


def set_precision(dps: int = 1200) -> None:
    mp.mp.dps = dps


def eye(n: int = 3) -> mp.matrix:
    return mp.eye(n)


def fro_norm(A: mp.matrix) -> mp.mpf:
    return mp.sqrt(mp.fsum(A[i, j] ** 2 for i in range(A.rows) for j in range(A.cols)))


def vnorm(v) -> mp.mpf:
    return mp.sqrt(mp.fsum(mp.mpf(x) ** 2 for x in v))


def dot(u, v) -> mp.mpf:
    return mp.fsum(mp.mpf(x) * mp.mpf(y) for x, y in zip(u, v))


def outer(u, v) -> mp.matrix:
    out = mp.zeros(len(u), len(v))
    for i, x in enumerate(u):
        for j, y in enumerate(v):
            out[i, j] = mp.mpf(x) * mp.mpf(y)
    return out


def cross(u, v):
    return [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]


def hat(v) -> mp.matrix:
    x, y, z = [mp.mpf(a) for a in v]
    return mp.matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])


def vee(S: mp.matrix):
    return [S[2, 1], S[0, 2], S[1, 0]]


def _sinc(theta: mp.mpf) -> mp.mpf:
    if abs(theta) < mp.mpf("1e-50"):
        t2 = theta * theta
        return mp.fsum([(-1) ** n * t2**n / mp.factorial(2 * n + 1) for n in range(30)])
    return mp.sin(theta) / theta


def _omc_over_t2(theta: mp.mpf) -> mp.mpf:
    if abs(theta) < mp.mpf("1e-50"):
        t2 = theta * theta
        return mp.fsum([(-1) ** n * t2**n / mp.factorial(2 * n + 2) for n in range(30)])
    return (1 - mp.cos(theta)) / (theta * theta)


def so3_exp_vec(omega) -> mp.matrix:
    omega = [mp.mpf(x) for x in omega]
    theta = vnorm(omega)
    K = hat(omega)
    return eye(3) + _sinc(theta) * K + _omc_over_t2(theta) * (K * K)


def trace3(A: mp.matrix) -> mp.mpf:
    return A[0, 0] + A[1, 1] + A[2, 2]


def so3_log_vec(R: mp.matrix):
    # atan2 form remains accurate close to the identity.
    S = mp.mpf("0.5") * (R - R.T)
    svec = vee(S)
    st = vnorm(svec)
    ct = (trace3(R) - 1) / 2
    ct = min(mp.mpf(1), max(mp.mpf(-1), ct))
    theta = mp.atan2(st, ct)
    if st == 0:
        return [mp.mpf(0), mp.mpf(0), mp.mpf(0)]
    return [theta * x / st for x in svec]


def det3(A: mp.matrix) -> mp.mpf:
    return (
        A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
        - A[0, 1] * (A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0])
        + A[0, 2] * (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0])
    )


def central_difference_operator(F_vec, M: mp.matrix, h: mp.mpf) -> mp.matrix:
    """Central-difference approximation D_{M,h} in the canonical basis of R^3."""
    if h <= 0:
        raise ValueError("h must be positive")
    D = mp.zeros(3, 3)
    for j in range(3):
        e = [mp.mpf(0)] * 3
        e[j] = mp.mpf(1)
        plus = M * so3_exp_vec([h * x for x in e])
        minus = M * so3_exp_vec([-h * x for x in e])
        fp = F_vec(plus)
        fm = F_vec(minus)
        for i in range(3):
            D[i, j] = (fp[i] - fm[i]) / (2 * h)
    return D


def _adapted_orthogonal_basis(a) -> mp.matrix:
    """Q=[q1 q2 q3], q1 parallel to a; columns are orthonormal."""
    theta = vnorm(a)
    if theta == 0:
        return eye(3)
    q1 = [x / theta for x in a]
    # Choose the coordinate direction least aligned with q1.
    idx = min(range(3), key=lambda i: abs(q1[i]))
    ref = [mp.mpf(0)] * 3
    ref[idx] = mp.mpf(1)
    q2 = cross(q1, ref)
    n2 = vnorm(q2)
    q2 = [x / n2 for x in q2]
    q3 = cross(q1, q2)
    n3 = vnorm(q3)
    q3 = [x / n3 for x in q3]
    Q = mp.zeros(3, 3)
    for i in range(3):
        Q[i, 0] = q1[i]
        Q[i, 1] = q2[i]
        Q[i, 2] = q3[i]
    return Q


def project_block_diagonal_so3(D: mp.matrix, a) -> mp.matrix:
    """Frobenius projection onto span{I, hat(a), hat(a)^2}.

    In a frame with first axis parallel to a, the projected matrix has the
    real 1+2 block form [[lambda,0,0],[0,alpha,-beta],[0,beta,alpha]].
    """
    if vnorm(a) == 0:
        return (trace3(D) / 3) * eye(3)
    Q = _adapted_orthogonal_basis(a)
    Dt = Q.T * D * Q
    lam = Dt[0, 0]
    alpha = (Dt[1, 1] + Dt[2, 2]) / 2
    beta = (Dt[2, 1] - Dt[1, 2]) / 2
    Bt = mp.matrix([[lam, 0, 0], [0, alpha, -beta], [0, beta, alpha]])
    return Q * Bt * Q.T


def rank_one_secant_operator(Dbase: mp.matrix, v, fX, fY) -> tuple[mp.matrix, list[mp.mpf]]:
    vv = dot(v, v)
    if vv == 0:
        raise ValueError("The two secant points must be distinct")
    vm = mp.matrix(v)
    delta = mp.matrix([fY[i] - fX[i] for i in range(3)])
    r = delta - Dbase * vm
    S = Dbase + outer([r[i] for i in range(3)], v) / vv
    return S, [r[i] for i in range(3)]


def bd_rank_one_divdiff(F_vec, X: mp.matrix, Y: mp.matrix, tau: mp.mpf = mp.mpf("0.1")):
    """Block-diagonal central-difference core + exact rank-one secant correction."""
    v = so3_log_vec(X.T * Y)
    nv = vnorm(v)
    if nv == 0:
        raise ValueError("The two secant points must be distinct")
    M = X * so3_exp_vec([mp.mpf("0.5") * x for x in v])
    h = tau * nv
    D = central_difference_operator(F_vec, M, h)
    Dbd = project_block_diagonal_so3(D, F_vec(M))
    S, r = rank_one_secant_operator(Dbd, v, F_vec(X), F_vec(Y))
    return S, v, h, D, Dbd, r


def r1_rank_one_divdiff(F_vec, X: mp.matrix, Y: mp.matrix, tau: mp.mpf = mp.mpf("0.1")):
    """General central-difference core + exact rank-one secant correction."""
    v = so3_log_vec(X.T * Y)
    nv = vnorm(v)
    if nv == 0:
        raise ValueError("The two secant points must be distinct")
    M = X * so3_exp_vec([mp.mpf("0.5") * x for x in v])
    h = tau * nv
    D = central_difference_operator(F_vec, M, h)
    S, r = rank_one_secant_operator(D, v, F_vec(X), F_vec(Y))
    return S, v, h, D, r


def computational_orders(residuals):
    rho = [mp.nan] * len(residuals)
    for k in range(1, len(residuals) - 1):
        den = mp.log(residuals[k] / residuals[k - 1])
        if den != 0:
            rho[k] = mp.log(residuals[k + 1] / residuals[k]) / den
    return rho


def sci(x, digits: int = 17) -> str:
    x = mp.mpf(x)
    if mp.isnan(x):
        return "--"
    if x == 0:
        return "0"
    return mp.nstr(x, digits, min_fixed=0, max_fixed=0)
