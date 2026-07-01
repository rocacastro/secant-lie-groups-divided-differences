"""
Restarted SO(3) secant example with high precision.

The new initial matrices X0 and X1 are the fifth and sixth iterates from the
previous corrected SO(3) experiment.  The computation uses mpmath with high
precision so that the table can display at least ten iterates before numerical
round-off dominates.

Run with:
    python so3_secant_example2_restart_high_precision.py
"""

import mpmath as mp
from _repo_utils import environment_text, results_dir, write_csv, write_text

# Working precision.  Increase this if you want to push the table further.
mp.mp.dps = 800

PAIRS = [(0, 1), (0, 2), (1, 2)]
TOL = mp.mpf("1e-700")
MIN_ROWS = 12       # print at least k = 0,...,11
MAX_INDEX = 12      # safety cap: compute no further than this index


def trace3(A):
    return sum(A[i, i] for i in range(3))


def fro_norm(A):
    return mp.sqrt(sum(A[i, j] ** 2 for i in range(3) for j in range(3)))


def vec_norm(v):
    return mp.sqrt(sum(v[i] ** 2 for i in range(3)))


def hat(v):
    """R^3 -> so(3)."""
    return mp.matrix(
        [
            [mp.mpf("0"), -v[2], v[1]],
            [v[2], mp.mpf("0"), -v[0]],
            [-v[1], v[0], mp.mpf("0")],
        ]
    )


def vee(S):
    """so(3) -> R^3."""
    return [S[2, 1], S[0, 2], S[1, 0]]


def sinc(theta):
    """sin(theta)/theta, evaluated by series near zero."""
    theta = mp.mpf(theta)
    if abs(theta) < mp.mpf("1e-100"):
        t2 = theta * theta
        s = mp.mpf("0")
        for n in range(40):
            s += (-1) ** n * t2**n / mp.factorial(2 * n + 1)
        return s
    return mp.sin(theta) / theta


def one_minus_cos_over_theta2(theta):
    """(1-cos(theta))/theta^2, evaluated by series near zero."""
    theta = mp.mpf(theta)
    if abs(theta) < mp.mpf("1e-100"):
        t2 = theta * theta
        s = mp.mpf("0")
        for n in range(40):
            s += (-1) ** n * t2**n / mp.factorial(2 * n + 2)
        return s
    return (1 - mp.cos(theta)) / (theta * theta)


def so3_exp_vec(omega):
    """Exponential map on SO(3), using Rodrigues' formula."""
    omega = [mp.mpf(x) for x in omega]
    theta = vec_norm(omega)
    K = hat(omega)
    return mp.eye(3) + sinc(theta) * K + one_minus_cos_over_theta2(theta) * (K * K)


def so3_log_vec(R):
    """Principal logarithm on SO(3), returned as a vector in R^3."""
    skew = mp.mpf("0.5") * (R - R.T)
    svec = vee(skew)
    sin_theta = vec_norm(svec)
    cos_theta = (trace3(R) - 1) / 2

    # Guard tiny numerical overshoots.
    if cos_theta > 1:
        cos_theta = mp.mpf("1")
    if cos_theta < -1:
        cos_theta = mp.mpf("-1")

    theta = mp.atan2(sin_theta, cos_theta)
    if sin_theta == 0:
        return [mp.mpf("0"), mp.mpf("0"), mp.mpf("0")]

    # Since sin_theta = ||vee((R-R^T)/2)||, this recovers omega robustly.
    coef = theta / sin_theta
    return [coef * x for x in svec]


def so3_log(R):
    return hat(so3_log_vec(R))


def F(X, Rtarget):
    """Vector field F(X)=log(X^T Rtarget)."""
    return so3_log(X.T * Rtarget)


def secant_so3(Rtarget, X0, X1, tol=TOL, min_rows=MIN_ROWS, max_index=MAX_INDEX):
    """
    Component-wise secant method on SO(3).

    The inverse divided difference is applied entrywise:
        W_ij = -F_k,ij * V_k,ij / (F_k,ij - F_{k-1},ij), i<j,
    where V_k = log(X_{k-1}^T X_k).
    """
    Xs = [X0.copy(), X1.copy()]
    Fs = [F(Xs[0], Rtarget), F(Xs[1], Rtarget)]
    residuals = [fro_norm(Fs[0]), fro_norm(Fs[1])]

    for k in range(1, max_index):
        V = so3_log(Xs[k - 1].T * Xs[k])
        Fprev = Fs[k - 1]
        Fcur = Fs[k]
        W = mp.zeros(3, 3)

        for i, j in PAIRS:
            denom = Fcur[i, j] - Fprev[i, j]
            if denom == 0:
                raise ZeroDivisionError(f"Zero divided-difference denominator at k={k}, entry {(i+1, j+1)}")
            W[i, j] = -Fcur[i, j] * V[i, j] / denom
            W[j, i] = -W[i, j]

        Xnext = Xs[k] * so3_exp_vec(vee(W))
        Xs.append(Xnext)
        Fs.append(F(Xnext, Rtarget))
        residuals.append(fro_norm(Fs[-1]))

        if len(residuals) >= min_rows and residuals[-1] < tol:
            break

    return Xs, residuals


def print_matrix(name, A, digits=13):
    print(f"{name} =")
    for i in range(3):
        row = "  ".join(mp.nstr(A[i, j], digits) for j in range(3))
        print("[" + row + "]")


def matrix_to_text(name, A, digits=20):
    lines = [f"{name} ="]
    for i in range(3):
        row = "  ".join(mp.nstr(A[i, j], digits) for j in range(3))
        lines.append("[" + row + "]")
    return "\n".join(lines)


def main():
    # Target rotation: angle pi/4 around the axis (1,1,1)/sqrt(3).
    axis = [1 / mp.sqrt(3), 1 / mp.sqrt(3), 1 / mp.sqrt(3)]
    theta = mp.pi / 4
    Rtarget = so3_exp_vec([theta * a for a in axis])

    # Previous corrected initial data used only to reproduce the old iterates.
    a0 = [mp.mpf("0.584379"), mp.mpf("-0.301106"), mp.mpf("0.440570")]
    a1 = [mp.mpf("0.589899"), mp.mpf("0.385493"), mp.mpf("0.291391")]
    Xold0 = Rtarget * so3_exp_vec(a0)
    Xold1 = Rtarget * so3_exp_vec(a1)

    # Compute old iterates up to X5 and X6, then restart from them.
    old_Xs, _ = secant_so3(Rtarget, Xold0, Xold1, min_rows=7, max_index=6)
    X0 = old_Xs[5]
    X1 = old_Xs[6]

    Xs, residuals = secant_so3(Rtarget, X0, X1)

    out = results_dir(__file__)
    write_csv(
        out / "so3_alignment_residuals.csv",
        ["k", "residual_norm_F"],
        [[k, mp.nstr(r, 40)] for k, r in enumerate(residuals)],
    )

    matrix_lines = [
        matrix_to_text("Rtarget", Rtarget),
        "",
        matrix_to_text("X0", X0),
        "det(X0) = " + mp.nstr(mp.det(X0), 40),
        "orthogonality_defect_X0 = " + mp.nstr(fro_norm(X0.T * X0 - mp.eye(3)), 40),
        "",
        matrix_to_text("X1", X1),
        "det(X1) = " + mp.nstr(mp.det(X1), 40),
        "orthogonality_defect_X1 = " + mp.nstr(fro_norm(X1.T * X1 - mp.eye(3)), 40),
        "",
        f"mp.dps = {mp.mp.dps}",
        f"tolerance = {mp.nstr(TOL, 10)}",
        environment_text(),
    ]
    write_text(out / "so3_alignment_input_matrices_and_summary.txt", "\n".join(matrix_lines) + "\n")

    output_lines = [
        matrix_to_text("Rtarget", Rtarget, 13),
        "",
        matrix_to_text("new X0 = old X5", X0, 13),
        "det(X0) = " + mp.nstr(mp.det(X0), 30),
        "||X0^T X0-I||_F = " + mp.nstr(fro_norm(X0.T * X0 - mp.eye(3)), 10),
        "",
        matrix_to_text("new X1 = old X6", X1, 13),
        "det(X1) = " + mp.nstr(mp.det(X1), 30),
        "||X1^T X1-I||_F = " + mp.nstr(fro_norm(X1.T * X1 - mp.eye(3)), 10),
        "",
        f"Tolerance used: epsilon = {mp.nstr(TOL, 5)}",
        "k    ||F(X_k)||_F",
    ]
    output_lines += [f"{k:2d}   {mp.nstr(r, 16)}" for k, r in enumerate(residuals)]
    output_lines += ["", f"Saved results to {out}"]
    write_text(out / "so3_alignment_output.txt", "\n".join(output_lines) + "\n")
    print("\n".join(output_lines))


if __name__ == "__main__":
    main()
