"""
Geodesic rotation averaging on SO(3) using the component-wise secant
Algorithm 3 from the manuscript.

The nonlinear equation is
    F(R) = sum_i w_i log(R^T R_i) = 0,
where F: SO(3) -> so(3). The divided difference is the component-wise
operator used in Algorithm 3:
    W_ij = - f_ij(R_k) (V_k)_ij / ( f_ij(R_k)-f_ij(R_{k-1}) ).

The code uses mpmath high precision and tolerance 1e-700.
"""

from __future__ import annotations
import mpmath as mp
from _repo_utils import environment_text, results_dir, write_csv, write_text

mp.mp.dps = 1200


def mat(rows):
    return mp.matrix([[mp.mpf(x) for x in row] for row in rows])


def eye(n=3):
    return mp.eye(n)


def fro_norm(A):
    return mp.sqrt(mp.fsum([A[i, j] ** 2 for i in range(A.rows) for j in range(A.cols)]))


def hat(v):
    x, y, z = [mp.mpf(a) for a in v]
    return mat([[0, -z, y], [z, 0, -x], [-y, x, 0]])


def vee(S):
    return [S[2, 1], S[0, 2], S[1, 0]]


def so3_exp_vec(omega):
    omega = [mp.mpf(x) for x in omega]
    theta = mp.sqrt(mp.fsum([x * x for x in omega]))
    K = hat(omega)
    I = eye(3)
    K2 = K * K
    if theta < mp.mpf('1e-80'):
        return I + K + mp.mpf('0.5') * K2 + (mp.mpf(1) / 6) * (K2 * K)
    return I + (mp.sin(theta) / theta) * K + ((1 - mp.cos(theta)) / (theta ** 2)) * K2


def so3_log_vec(R):
    tr = R[0, 0] + R[1, 1] + R[2, 2]
    c = (tr - 1) / 2
    c = max(min(c, mp.mpf(1)), mp.mpf(-1))
    theta = mp.acos(c)
    if theta < mp.mpf('1e-60'):
        return vee(mp.mpf('0.5') * (R - R.T))
    if mp.pi - theta < mp.mpf('1e-50'):
        raise ValueError('Principal logarithm is ill-conditioned near angle pi')
    S = (theta / (2 * mp.sin(theta))) * (R - R.T)
    return vee(S)


def so3_log(R):
    return hat(so3_log_vec(R))


def det3(A):
    return (A[0, 0] * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1])
            - A[0, 1] * (A[1, 0] * A[2, 2] - A[1, 2] * A[2, 0])
            + A[0, 2] * (A[1, 0] * A[2, 1] - A[1, 1] * A[2, 0]))


def F_rotation_average(R, rotations, weights):
    S = mp.zeros(3, 3)
    for w, Ri in zip(weights, rotations):
        S += w * so3_log(R.T * Ri)
    return S


def objective(R, rotations, weights):
    return mp.mpf('0.5') * mp.fsum([
        w * fro_norm(so3_log(R.T * Ri)) ** 2 for w, Ri in zip(weights, rotations)
    ])


def secant_algorithm3_componentwise(R0, R1, rotations, weights, tol=mp.mpf('1e-700'), max_iter=50):
    Xs = [mp.matrix(R0), mp.matrix(R1)]
    Fs = [F_rotation_average(Xs[0], rotations, weights),
          F_rotation_average(Xs[1], rotations, weights)]

    for k in range(1, max_iter):
        Xprev, Xcur = Xs[k - 1], Xs[k]
        Fprev, Fcur = Fs[k - 1], Fs[k]

        # X_cur = X_prev exp(V_k)
        V = so3_log(Xprev.T * Xcur)
        W = mp.zeros(3, 3)

        # Component-wise divided difference used in Algorithm 3.
        for i, j in [(0, 1), (0, 2), (1, 2)]:
            denom = Fcur[i, j] - Fprev[i, j]
            if abs(denom) < mp.mpf('1e-1100'):
                raise ZeroDivisionError(f'Small denominator at k={k}, entry ({i+1},{j+1})')
            W[i, j] = -Fcur[i, j] * V[i, j] / denom
            W[j, i] = -W[i, j]

        Xnext = Xcur * so3_exp_vec(vee(W))
        Fnext = F_rotation_average(Xnext, rotations, weights)
        Xs.append(Xnext)
        Fs.append(Fnext)
        if fro_norm(Fnext) < tol:
            break

    return Xs, Fs


def sci(x, digits=17):
    x = mp.mpf(x)
    if x == 0:
        return '0'
    return mp.nstr(x, digits, min_fixed=0, max_fixed=0)


def matrix_to_text(name, A, digits=20):
    lines = [f"{name} ="]
    for i in range(3):
        row = "  ".join(mp.nstr(A[i, j], digits) for j in range(3))
        lines.append("[" + row + "]")
    return "\n".join(lines)


def main():
    # Nominal unknown mean rotation.
    r_star = [mp.mpf('0.45'), mp.mpf('-0.25'), mp.mpf('0.35')]
    R_star = so3_exp_vec(r_star)

    # Deterministic noisy rotation measurements around R_star.
    # The perturbations model small roll/pitch/yaw measurement errors.
    perturbations = []
    for a, e in [
        (mp.mpf('0.060'), [1, 0, 0]),
        (mp.mpf('0.050'), [0, 1, 0]),
        (mp.mpf('0.040'), [0, 0, 1]),
    ]:
        perturbations.append([a * ee for ee in e])
        perturbations.append([-a * ee for ee in e])

    weights = [mp.mpf(1) / len(perturbations)] * len(perturbations)
    rotations = [R_star * so3_exp_vec(xi) for xi in perturbations]

    # Two initial rotations in SO(3). They are close enough to display the local regime.
    R0 = R_star * so3_exp_vec([mp.mpf('0.050'), mp.mpf('-0.040'), mp.mpf('0.030')])
    R1 = R_star * so3_exp_vec([mp.mpf('-0.040'), mp.mpf('0.030'), mp.mpf('0.020')])

    Xs, Fs = secant_algorithm3_componentwise(R0, R1, rotations, weights)
    out = results_dir(__file__)

    table_rows = []
    for k, (X, FX) in enumerate(zip(Xs, Fs)):
        table_rows.append([
            k,
            sci(fro_norm(FX), 40),
            sci(objective(X, rotations, weights), 40),
            sci(fro_norm(X.T * X - eye(3)), 40),
            sci(det3(X), 40),
        ])
    write_csv(
        out / "so3_rotation_averaging_table.csv",
        ["k", "residual_norm_F", "objective_J", "orthogonality_defect", "determinant"],
        table_rows,
    )

    input_lines = [
        f"mp.dps = {mp.mp.dps}",
        "tolerance = 1e-700",
        "r_star = (0.45, -0.25, 0.35)",
        "weights = " + ", ".join(mp.nstr(w, 30) for w in weights),
        "perturbations:",
    ]
    for i, xi in enumerate(perturbations, start=1):
        input_lines.append(f"xi_{i} = (" + ", ".join(mp.nstr(x, 30) for x in xi) + ")")
    input_lines += [
        "",
        matrix_to_text("R_star", R_star),
        "",
        matrix_to_text("R0", R0),
        "",
        matrix_to_text("R1", R1),
        "",
        "Generated rotations R_i:",
    ]
    for i, Ri in enumerate(rotations, start=1):
        input_lines += ["", matrix_to_text(f"R_{i}", Ri)]
    input_lines += ["", "Environment", environment_text()]
    write_text(out / "so3_rotation_averaging_inputs.txt", "\n".join(input_lines) + "\n")

    lines = [
        f"mp.dps = {mp.mp.dps}",
        "tolerance = 1e-700",
        "J(R_star) = " + sci(objective(R_star, rotations, weights), 20),
        "||F(R_star)||_F = " + sci(fro_norm(F_rotation_average(R_star, rotations, weights)), 20),
        "det(R_star) = " + sci(det3(R_star), 20),
        "",
        "Table",
        "k  ||F(R_k)||_F             J(R_k)                  ||R_k^T R_k-I||_F",
    ]
    for k, res, obj, ortho, _det in table_rows:
        lines.append(f"{int(k):2d}  {res:>22}  {obj:>22}  {ortho:>22}")
    lines += [
        "",
        "final residual = " + sci(fro_norm(Fs[-1]), 30),
        "iterations = " + str(len(Xs) - 1),
        "det(R_final) = " + sci(det3(Xs[-1]), 30),
        "",
        f"Saved results to {out}",
    ]
    write_text(out / "so3_rotation_averaging_output.txt", "\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == '__main__':
    main()
