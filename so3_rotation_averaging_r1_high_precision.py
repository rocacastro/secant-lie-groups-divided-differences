"""Example 3: geodesic rotation averaging using the general R1 divided difference.

The auxiliary operator is the central-difference approximation at the intrinsic
midpoint; a rank-one correction imposes the secant equation exactly.  No
projection is used because the derivative of this field is not generally in the
SO(3) alignment subspace span{I,hat(f(M)),hat(f(M))^2}.
"""
from __future__ import annotations
import mpmath as mp
from _repo_utils import environment_text
from so3_rankone_divdiff_utils import (
    set_precision, results_dir, write_csv, write_text, eye, fro_norm, vnorm,
    hat, so3_exp_vec, so3_log_vec, r1_rank_one_divdiff, computational_orders,
    sci, det3,
)

DPS = 1200
set_precision(DPS)
TOL = mp.mpf("1e-700")
MAX_ITER = 40
TAU = mp.mpf("0.1")


def F_rotation_average_factory(rotations, weights):
    def F_vec(R):
        out = [mp.mpf(0), mp.mpf(0), mp.mpf(0)]
        for w, Ri in zip(weights, rotations):
            q = so3_log_vec(R.T * Ri)
            for j in range(3):
                out[j] += w * q[j]
        return out
    return F_vec


def objective(R, rotations, weights):
    return mp.mpf("0.5") * mp.fsum(
        w * fro_norm(hat(so3_log_vec(R.T * Ri))) ** 2
        for w, Ri in zip(weights, rotations)
    )


def solve(R0, R1, rotations, weights):
    F_vec = F_rotation_average_factory(rotations, weights)
    Xs = [mp.matrix(R0), mp.matrix(R1)]
    fs = [F_vec(R0), F_vec(R1)]
    res = [fro_norm(hat(fs[0])), fro_norm(hat(fs[1]))]
    secant_defects = [mp.nan, mp.nan]
    inverse_norms = [mp.nan, mp.nan]

    for k in range(1, MAX_ITER):
        S, v, h, D, r = r1_rank_one_divdiff(F_vec, Xs[k - 1], Xs[k], TAU)
        stepm = mp.lu_solve(S, -mp.matrix(fs[k]))
        step = [stepm[i] for i in range(3)]
        delta = mp.matrix([fs[k][i] - fs[k - 1][i] for i in range(3)])
        secant_defects.append(vnorm([x for x in (S * mp.matrix(v) - delta)]))
        inverse_norms.append(fro_norm(S ** -1))

        Xnext = Xs[k] * so3_exp_vec(step)
        Xs.append(Xnext)
        fs.append(F_vec(Xnext))
        res.append(fro_norm(hat(fs[-1])))
        if res[-1] < TOL:
            break
    return Xs, fs, res, secant_defects, inverse_norms


def main():
    r_star = [mp.mpf("0.45"), mp.mpf("-0.25"), mp.mpf("0.35")]
    R_star = so3_exp_vec(r_star)
    perturbations = []
    for amp, e in [
        (mp.mpf("0.060"), [1, 0, 0]),
        (mp.mpf("0.050"), [0, 1, 0]),
        (mp.mpf("0.040"), [0, 0, 1]),
    ]:
        perturbations.append([amp * x for x in e])
        perturbations.append([-amp * x for x in e])
    weights = [mp.mpf(1) / 6] * 6
    rotations = [R_star * so3_exp_vec(xi) for xi in perturbations]

    R0 = R_star * so3_exp_vec([mp.mpf("0.050"), mp.mpf("-0.040"), mp.mpf("0.030")])
    R1 = R_star * so3_exp_vec([mp.mpf("-0.040"), mp.mpf("0.030"), mp.mpf("0.020")])

    Xs, fs, res, defects, invnorms = solve(R0, R1, rotations, weights)
    rho = computational_orders(res)
    out = results_dir(__file__)

    rows = []
    for k in range(len(res)):
        rows.append([
            k, sci(res[k], 42), sci(objective(Xs[k], rotations, weights), 42),
            sci(rho[k], 24), sci(fro_norm(Xs[k].T * Xs[k] - eye(3)), 30),
            sci(det3(Xs[k]), 30), sci(defects[k], 20), sci(invnorms[k], 20),
        ])
    write_csv(
        out / "so3_rotation_averaging_r1_table.csv",
        ["k", "residual_norm_F", "objective_J", "rho_k", "orthogonality_defect",
         "determinant", "secant_defect", "inverse_frobenius_norm"],
        rows,
    )

    lines = [
        "SO(3) geodesic rotation averaging with central-difference core + rank-one correction",
        f"mp.dps = {mp.mp.dps}", f"tau = {TAU}", f"tolerance = {TOL}",
        "k  ||F(R_k)||_F                  J(R_k)                    rho_k",
    ]
    for k in range(len(res)):
        lines.append(f"{k:2d}  {sci(res[k], 22):>26}  {sci(objective(Xs[k], rotations, weights), 22):>26}  {sci(rho[k], 12):>14}")
    lines += [
        "", "final residual = " + sci(res[-1], 40),
        "secant updates = " + str(len(res) - 2),
        "final orthogonality defect = " + sci(fro_norm(Xs[-1].T * Xs[-1] - eye(3)), 30),
        "max secant defect = " + sci(max([d for d in defects[2:] if not mp.isnan(d)]), 30),
        "", "Environment", environment_text(),
        "saved to " + str(out),
    ]
    write_text(out / "so3_rotation_averaging_r1_output.txt", "\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
