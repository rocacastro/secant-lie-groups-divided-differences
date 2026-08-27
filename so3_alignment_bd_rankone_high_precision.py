"""Example 2: SO(3) orientation alignment using the new BD divided difference.

The auxiliary operator is obtained by central differences at the intrinsic
midpoint, projected onto span{I, hat(f(M)), hat(f(M))^2}, and corrected by an
exact rank-one secant term.  The algorithm evaluates no derivative of F.
"""
from __future__ import annotations
import mpmath as mp
from _repo_utils import environment_text
from so3_rankone_divdiff_utils import (
    set_precision, results_dir, write_csv, write_text, eye, fro_norm, vnorm,
    hat, so3_exp_vec, so3_log_vec, bd_rank_one_divdiff, computational_orders,
    sci, det3,
)

DPS = 1200
set_precision(DPS)
TOL = mp.mpf("1e-700")
MAX_ITER = 30
TAU = mp.mpf("0.1")


def F_factory(Rtarget):
    def F_vec(X):
        return so3_log_vec(X.T * Rtarget)
    return F_vec


def solve(Rtarget, X0, X1):
    F_vec = F_factory(Rtarget)
    Xs = [mp.matrix(X0), mp.matrix(X1)]
    fs = [F_vec(X0), F_vec(X1)]
    res = [fro_norm(hat(fs[0])), fro_norm(hat(fs[1]))]
    secant_defects = [mp.nan, mp.nan]
    inverse_norms = [mp.nan, mp.nan]
    correction_norms = [mp.nan, mp.nan]

    for k in range(1, MAX_ITER):
        S, v, h, D, Dbd, r = bd_rank_one_divdiff(F_vec, Xs[k - 1], Xs[k], TAU)
        stepm = mp.lu_solve(S, -mp.matrix(fs[k]))
        step = [stepm[i] for i in range(3)]

        delta = mp.matrix([fs[k][i] - fs[k - 1][i] for i in range(3)])
        defect = vnorm([x for x in (S * mp.matrix(v) - delta)])
        secant_defects.append(defect)
        inverse_norms.append(fro_norm(S ** -1))
        correction_norms.append(vnorm(r) / vnorm(v))

        Xnext = Xs[k] * so3_exp_vec(step)
        Xs.append(Xnext)
        fs.append(F_vec(Xnext))
        res.append(fro_norm(hat(fs[-1])))
        if res[-1] < TOL:
            break

    return Xs, fs, res, secant_defects, inverse_norms, correction_norms


def main():
    axis = [1 / mp.sqrt(3)] * 3
    theta = mp.pi / 4
    Rtarget = so3_exp_vec([theta * x for x in axis])

    b0 = [mp.mpf("-0.1309040431"), mp.mpf("-0.3172408338"), mp.mpf("0.0427275100")]
    b1 = [mp.mpf("-0.3234918357"), mp.mpf("0.0015621005"), mp.mpf("-0.0459871328")]
    X0 = Rtarget * so3_exp_vec(b0)
    X1 = Rtarget * so3_exp_vec(b1)

    Xs, fs, res, defects, invnorms, corrnorms = solve(Rtarget, X0, X1)
    rho = computational_orders(res)
    out = results_dir(__file__)

    rows = []
    for k in range(len(res)):
        rows.append([
            k, sci(res[k], 42), sci(rho[k], 24),
            sci(fro_norm(Xs[k].T * Xs[k] - eye(3)), 30),
            sci(det3(Xs[k]), 30), sci(defects[k], 20),
            sci(invnorms[k], 20), sci(corrnorms[k], 20),
        ])
    write_csv(
        out / "so3_alignment_bd_rankone_table.csv",
        ["k", "residual_norm_F", "rho_k", "orthogonality_defect", "determinant",
         "secant_defect", "inverse_frobenius_norm", "rank_one_correction_proxy"],
        rows,
    )

    lines = [
        "SO(3) alignment with block-diagonal central-difference core + rank-one correction",
        f"mp.dps = {mp.mp.dps}", f"tau = {TAU}", f"tolerance = {TOL}",
        "k  ||F(X_k)||_F                  rho_k          orthogonality defect",
    ]
    for k in range(len(res)):
        lines.append(f"{k:2d}  {sci(res[k], 22):>26}  {sci(rho[k], 12):>14}  {sci(fro_norm(Xs[k].T*Xs[k]-eye(3)), 12):>20}")
    lines += [
        "", "final residual = " + sci(res[-1], 40),
        "secant updates = " + str(len(res) - 2),
        "final orthogonality defect = " + sci(fro_norm(Xs[-1].T * Xs[-1] - eye(3)), 30),
        "max secant defect = " + sci(max([d for d in defects[2:] if not mp.isnan(d)]), 30),
        "", "Environment", environment_text(),
        "saved to " + str(out),
    ]
    write_text(out / "so3_alignment_bd_rankone_output.txt", "\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
