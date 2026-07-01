"""
Benchmark comparing Algorithm 1 and Algorithm 2 on SO(2).

This script reproduces the residual comparison and the CPU-time comparison used
in the manuscript. It saves CSV/TXT output in results/current/ next to this
script. Timing values depend on hardware, Python version and system load; the
residuals and order estimates are deterministic.
"""
from __future__ import annotations

import mpmath as mp
from time import perf_counter, process_time
from _repo_utils import environment_text, results_dir, write_csv, write_text

mp.mp.dps = 2500
SQRT2 = mp.sqrt(2)
ONE_QUARTER = mp.mpf(1) / 4
MAX_ITER = 14
REPEATS = 100


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
    return [[mp.mpf(0), -1 / v], [1 / v, mp.mpf(0)]]


def f12(P):
    return P[0][1] - P[1][0] + ONE_QUARTER


def residual_from_entry(f):
    return SQRT2 * abs(f)


def algorithm1(max_iter=MAX_ITER):
    v0 = mp.mpf(1) / 7
    P_prev = [[mp.mpf(1), mp.mpf(0)], [mp.mpf(0), mp.mpf(1)]]
    P_cur = rot(v0)
    v_prev = v0
    residuals = [frobenius_norm(F(P_prev)), frobenius_norm(F(P_cur))]

    for _ in range(1, max_iter):
        G = matmul(matsub(P_cur, P_prev), inv_skew(v_prev))
        trace_G = G[0][0] + G[1][1]
        F_cur = F(P_cur)
        V_new = [[-F_cur[i][j] / trace_G for j in range(2)] for i in range(2)]
        v_new = V_new[0][1]
        P_next = matmul(P_cur, rot(v_new))
        residuals.append(frobenius_norm(F(P_next)))
        P_prev, P_cur, v_prev = P_cur, P_next, v_new

    return residuals


def algorithm2(max_iter=MAX_ITER):
    v0 = mp.mpf(1) / 7
    P_prev = [[mp.mpf(1), mp.mpf(0)], [mp.mpf(0), mp.mpf(1)]]
    P_cur = rot(v0)
    v_prev = v0
    b = f12(P_prev)
    c = f12(P_cur)
    residuals = [residual_from_entry(b), residual_from_entry(c)]

    for _ in range(1, max_iter):
        v_new = (c / (b - c)) * v_prev
        P_next = matmul(P_cur, rot(v_new))
        f_next = f12(P_next)
        residuals.append(residual_from_entry(f_next))
        P_prev, P_cur, v_prev = P_cur, P_next, v_new
        b, c = c, f_next

    return residuals


def convergence_order(residuals, k):
    return mp.log(residuals[k + 1] / residuals[k]) / mp.log(residuals[k] / residuals[k - 1])


def timed_average(func, repeats=REPEATS):
    func()  # warm-up
    t0 = perf_counter()
    c0 = process_time()
    for _ in range(repeats):
        func()
    wall = perf_counter() - t0
    cpu = process_time() - c0
    return wall / repeats, cpu / repeats


def main():
    r1 = algorithm1()
    r2 = algorithm2()
    out = results_dir(__file__)

    comparison_rows = [
        [k, mp.nstr(a, 30), mp.nstr(b, 30), mp.nstr(abs(a - b), 30)]
        for k, (a, b) in enumerate(zip(r1, r2))
    ]
    write_csv(out / "so2_residual_comparison_alg1_alg2.csv",
              ["k", "algorithm1_residual", "algorithm2_residual", "abs_difference"],
              comparison_rows)

    order_rows = []
    orders = []
    for k in range(9, 14):
        rho1 = convergence_order(r1, k)
        rho2 = convergence_order(r2, k)
        orders.append(rho1)
        order_rows.append([k, mp.nstr(rho1, 30), mp.nstr(rho2, 30)])
    mean_order = sum(orders) / len(orders)
    write_csv(out / "so2_order_estimates_alg1_alg2.csv",
              ["k", "rho_algorithm1", "rho_algorithm2"], order_rows)

    wall1, cpu1 = timed_average(algorithm1)
    wall2, cpu2 = timed_average(algorithm2)
    efficiency_rows = [
        ["Algorithm 1", MAX_ITER - 1, mp.nstr(mean_order, 20), f"{wall1:.10e}", f"{cpu1:.10e}", "1.0", "1.0"],
        ["Algorithm 2", MAX_ITER - 1, mp.nstr(mean_order, 20), f"{wall2:.10e}", f"{cpu2:.10e}", f"{wall1/wall2:.10e}", f"{cpu1/cpu2:.10e}"],
    ]
    write_csv(out / "so2_efficiency_benchmark.csv",
              ["method", "iterations", "mean_order", "wall_time_seconds", "cpu_time_seconds", "relative_wall_efficiency", "relative_cpu_efficiency"],
              efficiency_rows)

    lines = []
    lines += ["Residual comparison", "k  Algorithm 1              Algorithm 2              abs difference"]
    for k, a, b, diff in comparison_rows:
        lines.append(f"{int(k):2d} {a:>24} {b:>24} {diff:>16}")
    lines += ["", "Observed computational order"]
    for k, rho1, rho2 in order_rows:
        lines.append(f"k={int(k):2d}: rho1={mp.nstr(mp.mpf(rho1), 12)}, rho2={mp.nstr(mp.mpf(rho2), 12)}")
    lines += [
        "",
        f"Mean observed order = {mp.nstr(mean_order, 14)}",
        "",
        "Timing benchmark",
        f"mp.dps  = {mp.mp.dps}",
        f"repeats = {REPEATS}",
        f"Algorithm 1: wall/run = {wall1:.8f} s, cpu/run = {cpu1:.8f} s",
        f"Algorithm 2: wall/run = {wall2:.8f} s, cpu/run = {cpu2:.8f} s",
        f"Relative wall efficiency T1/T2 = {wall1/wall2:.4f}",
        f"Relative CPU efficiency  T1/T2 = {cpu1/cpu2:.4f}",
        "",
        "Environment",
        environment_text(),
    ]
    write_text(out / "so2_algorithms_efficiency_benchmark_output.txt", "\n".join(lines) + "\n")

    print("\n".join(lines))
    print(f"\nSaved results to {out}")


if __name__ == "__main__":
    main()
