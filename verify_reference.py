"""Verify deterministic numerical columns against the archived reference outputs.

Run ``python run_all.py`` first. Timing columns are intentionally not checked,
because they depend on hardware and system load.
"""
from __future__ import annotations

import csv
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parent
CURRENT = ROOT / "results" / "current"
REFERENCE = ROOT / "results" / "reference"
mp.mp.dps = 80


class VerificationError(RuntimeError):
    pass


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise VerificationError(f"Missing file: {path}")
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _number(text: str) -> mp.mpf | None:
    value = text.strip()
    if value in {"", "--", "nan", "NaN"}:
        return None
    return mp.mpf(value)


def _close(a: mp.mpf, b: mp.mpf, rtol: mp.mpf, atol: mp.mpf) -> bool:
    return abs(a - b) <= atol + rtol * max(abs(a), abs(b))


def compare_csv(
    generated_name: str,
    reference_name: str,
    columns: list[str],
    *,
    rtol: str = "1e-24",
    atol: str = "1e-80",
) -> None:
    generated = _read_csv(CURRENT / generated_name)
    reference = _read_csv(REFERENCE / reference_name)
    if len(generated) != len(reference):
        raise VerificationError(
            f"Row-count mismatch for {generated_name}: "
            f"{len(generated)} != {len(reference)}"
        )
    rtol_m = mp.mpf(rtol)
    atol_m = mp.mpf(atol)
    for row_index, (g_row, r_row) in enumerate(zip(generated, reference)):
        for column in columns:
            ga = _number(g_row[column])
            rb = _number(r_row[column])
            if ga is None or rb is None:
                if ga != rb:
                    raise VerificationError(
                        f"Missing-value mismatch in {generated_name}, "
                        f"row {row_index}, column {column}"
                    )
                continue
            if not _close(ga, rb, rtol_m, atol_m):
                raise VerificationError(
                    f"Mismatch in {generated_name}, row {row_index}, column {column}: "
                    f"{ga} != {rb}"
                )
    print(f"PASS: {generated_name}")


def main() -> None:
    compare_csv(
        "so2_algorithm1_residuals.csv",
        "so2_a1_res.csv",
        ["k", "residual_norm_F"],
    )
    compare_csv(
        "so2_algorithm1_orders.csv",
        "so2_a1_ord.csv",
        ["k", "rho_k"],
    )
    compare_csv(
        "so2_algorithm2_residuals.csv",
        "so2_a2_res.csv",
        ["k", "residual_norm_F"],
    )
    compare_csv(
        "so2_algorithm2_orders.csv",
        "so2_a2_ord.csv",
        ["k", "rho_k"],
    )
    compare_csv(
        "so2_residual_comparison_alg1_alg2.csv",
        "so2_cmp_res.csv",
        ["k", "algorithm1_residual", "algorithm2_residual"],
    )
    compare_csv(
        "so2_order_estimates_alg1_alg2.csv",
        "so2_cmp_ord.csv",
        ["k", "rho_algorithm1", "rho_algorithm2"],
    )
    compare_csv(
        "so3_alignment_bd_rankone_table.csv",
        "so3_align.csv",
        ["k", "residual_norm_F", "rho_k", "orthogonality_defect", "determinant"],
        rtol="1e-20",
        atol="1e-1190",
    )
    compare_csv(
        "so3_rotation_averaging_r1_table.csv",
        "so3_avg.csv",
        ["k", "residual_norm_F", "objective_J", "rho_k", "orthogonality_defect", "determinant"],
        rtol="1e-20",
        atol="1e-1190",
    )
    print("\nAll deterministic reference checks passed.")


if __name__ == "__main__":
    main()
