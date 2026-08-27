"""Run all reproducibility scripts for the numerical experiments in the manuscript."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    "so2_algorithm1_matrix_high_precision.py",
    "so2_algorithm2_optimized_high_precision.py",
    "so2_algorithms_efficiency_benchmark.py",
    "so3_alignment_bd_rankone_high_precision.py",
    "so3_rotation_averaging_r1_high_precision.py",
]


def main() -> None:
    for script in SCRIPTS:
        path = ROOT / script
        if not path.is_file():
            raise FileNotFoundError(f"Missing reproducibility script: {path}")
        print("=" * 80)
        print(f"Running {script}")
        print("=" * 80)
        subprocess.run([sys.executable, str(path)], check=True, cwd=ROOT)
    print("\nAll scripts completed successfully.")
    print(f"Generated results: {ROOT / 'results' / 'current'}")


if __name__ == "__main__":
    main()
