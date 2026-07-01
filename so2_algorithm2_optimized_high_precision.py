"""
Algorithm 2 on SO(2) with high-precision arithmetic.

This script reproduces the Algorithm 2 residuals and order estimates used in the
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
    """Return exp(vJ), J=[[0,1],[-1,0]]."""
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


def f12(P):
    """(1,2)-entry of F(P)=P-P^T+A, A12=1/4."""
    return P[0][1] - P[1][0] + ONE_QUARTER


def residual_from_entry(f):
    """Frobenius norm of [[0,f],[-f,0]]."""
    return SQRT2 * abs(f)


def algorithm2(max_iter=MAX_ITER):
    v0 = mp.mpf(1) / 7
    P_prev = [[mp.mpf(1), mp.mpf(0)], [mp.mpf(0), mp.mpf(1)]]
    P_cur = rot(v0)
    v_prev = v0

    b = f12(P_prev)
    c = f12(P_cur)
    residuals = [residual_from_entry(b), residual_from_entry(c)]
    iterates = [P_prev, P_cur]
    steps = [v_prev]

    for _ in range(1, max_iter):
        # Component-wise divided difference formula:
        # V_n=(c_12/(b_12-c_12)) V_{n-1}.
        v_new = (c / (b - c)) * v_prev
        P_next = matmul(P_cur, rot(v_new))
        f_next = f12(P_next)

        iterates.append(P_next)
        steps.append(v_new)
        residuals.append(residual_from_entry(f_next))
        P_prev, P_cur, v_prev = P_cur, P_next, v_new
        b, c = c, f_next

    return residuals, iterates, steps


def convergence_order(residuals, k):
    return mp.log(residuals[k + 1] / residuals[k]) / mp.log(residuals[k] / residuals[k - 1])


def main():
    residuals, _, _ = algorithm2()
    orders = [(k, convergence_order(residuals, k)) for k in range(9, 14)]
    out = results_dir(__file__)

    write_csv(
        out / "so2_algorithm2_residuals.csv",
        ["k", "residual_norm_F"],
        [[k, mp.nstr(r, 30)] for k, r in enumerate(residuals)],
    )
    write_csv(
        out / "so2_algorithm2_orders.csv",
        ["k", "rho_k"],
        [[k, mp.nstr(rho, 30)] for k, rho in orders],
    )

    lines = [
        "Algorithm 2 on SO(2)",
        f"mp.dps = {mp.mp.dps}",
        f"max_iter = {MAX_ITER}",
        environment_text(),
        "",
        "Residual norms:",
    ]
    lines += [f"{k:2d}  {mp.nstr(r, 16)}" for k, r in enumerate(residuals)]
    lines += ["", "Observed computational order:"]
    lines += [f"rho_{k} = {mp.nstr(rho, 12)}" for k, rho in orders]
    write_text(out / "so2_algorithm2_summary.txt", "\n".join(lines) + "\n")

    print("Algorithm 2: residual norms")
    for k, r in enumerate(residuals):
        print(f"{k:2d}  {mp.nstr(r, 16)}")
    print("\nAlgorithm 2: observed computational order")
    for k, rho in orders:
        print(f"rho_{k} = {mp.nstr(rho, 12)}")
    print(f"\nSaved results to {out}")


if __name__ == "__main__":
    main()
