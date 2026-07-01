"""
Algorithm 1 on SO(2) with high-precision arithmetic.

This script reproduces the Algorithm 1 residuals and order estimates used in the
SO(2) numerical comparison of the manuscript. It also saves CSV/TXT output in
results/current/ next to this script.
"""
from __future__ import annotations

import mpmath as mp
from _repo_utils import environment_text, results_dir, write_csv, write_text

mp.mp.dps = 2500
SQRT2 = mp.sqrt(2)
ONE_QUARTER = mp.mpf(1) / 4
MAX_ITER = 14


def rot(v):
    c = mp.cos(v)
    s = mp.sin(v)
    return [[c, s], [-s, c]]


def matmul(A, B):
    return [
        [A[0][0] * B[0][0] + A[0][1] * B[1][0],
         A[0][0] * B[0][1] + A[0][1] * B[1][1]],
        [A[1][0] * B[0][0] + A[1][1] * B[1][0],
         A[1][0] * B[0][1] + A[1][1] * B[1][1]],
    ]


def matsub(A, B):
    return [[A[0][0] - B[0][0], A[0][1] - B[0][1]],
            [A[1][0] - B[1][0], A[1][1] - B[1][1]]]


def transpose(A):
    return [[A[0][0], A[1][0]], [A[0][1], A[1][1]]]


def F(P):
    PT = transpose(P)
    return [[P[0][0] - PT[0][0], P[0][1] - PT[0][1] + ONE_QUARTER],
            [P[1][0] - PT[1][0] - ONE_QUARTER, P[1][1] - PT[1][1]]]


def frobenius_norm(M):
    return mp.sqrt(sum(M[i][j] ** 2 for i in range(2) for j in range(2)))


def inv_skew(v):
    """Inverse of V=vJ, J=[[0,1],[-1,0]]."""
    return [[mp.mpf(0), -1 / v], [1 / v, mp.mpf(0)]]


def algorithm1(max_iter=MAX_ITER):
    v0 = mp.mpf(1) / 7
    P_prev = [[mp.mpf(1), mp.mpf(0)], [mp.mpf(0), mp.mpf(1)]]
    P_cur = rot(v0)
    v_prev = v0

    residuals = [frobenius_norm(F(P_prev)), frobenius_norm(F(P_cur))]
    iterates = [P_prev, P_cur]
    steps = [v_prev]

    for _ in range(1, max_iter):
        # Optimized equivalent of P_{n-1}(exp(V_{n-1})-I)V_{n-1}^{-1}.
        G = matmul(matsub(P_cur, P_prev), inv_skew(v_prev))
        trace_G = G[0][0] + G[1][1]
        F_cur = F(P_cur)

        V_new = [[-F_cur[i][j] / trace_G for j in range(2)] for i in range(2)]
        v_new = V_new[0][1]
        P_next = matmul(P_cur, rot(v_new))

        iterates.append(P_next)
        steps.append(v_new)
        residuals.append(frobenius_norm(F(P_next)))
        P_prev, P_cur, v_prev = P_cur, P_next, v_new

    return residuals, iterates, steps


def convergence_order(residuals, k):
    return mp.log(residuals[k + 1] / residuals[k]) / mp.log(residuals[k] / residuals[k - 1])


def main():
    residuals, _, _ = algorithm1()
    orders = [(k, convergence_order(residuals, k)) for k in range(9, 14)]
    out = results_dir(__file__)

    write_csv(
        out / "so2_algorithm1_residuals.csv",
        ["k", "residual_norm_F"],
        [[k, mp.nstr(r, 30)] for k, r in enumerate(residuals)],
    )
    write_csv(
        out / "so2_algorithm1_orders.csv",
        ["k", "rho_k"],
        [[k, mp.nstr(rho, 30)] for k, rho in orders],
    )

    lines = [
        "Algorithm 1 on SO(2)",
        f"mp.dps = {mp.mp.dps}",
        f"max_iter = {MAX_ITER}",
        environment_text(),
        "",
        "Residual norms:",
    ]
    lines += [f"{k:2d}  {mp.nstr(r, 16)}" for k, r in enumerate(residuals)]
    lines += ["", "Observed computational order:"]
    lines += [f"rho_{k} = {mp.nstr(rho, 12)}" for k, rho in orders]
    write_text(out / "so2_algorithm1_summary.txt", "\n".join(lines) + "\n")

    print("Algorithm 1: residual norms")
    for k, r in enumerate(residuals):
        print(f"{k:2d}  {mp.nstr(r, 16)}")
    print("\nAlgorithm 1: observed computational order")
    for k, rho in orders:
        print(f"rho_{k} = {mp.nstr(rho, 12)}")
    print(f"\nSaved results to {out}")


if __name__ == "__main__":
    main()
