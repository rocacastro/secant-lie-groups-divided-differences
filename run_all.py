"""Run all reproducibility scripts for the numerical examples in the manuscript."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = [
    "so2_algorithm1_matrix_high_precision.py",
    "so2_algorithm2_optimized_high_precision.py",
    "so2_algorithms_efficiency_benchmark.py",
    "so3_secant_example2_restart_high_precision.py",
    "so3_rotation_averaging_componentwise_alg3_high_precision.py",
]


def main() -> None:
    for script in SCRIPTS:
        print("=" * 80)
        print(f"Running {script}")
        print("=" * 80)
        subprocess.run([sys.executable, str(ROOT / script)], check=True)
    print("\nAll scripts completed. Results are in results/current/.")


if __name__ == "__main__":
    main()
