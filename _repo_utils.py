"""Small I/O helpers for the reproducibility scripts."""
from __future__ import annotations

import csv
import platform
from pathlib import Path
from datetime import datetime, timezone


def package_root(current_file: str) -> Path:
    """Return the directory containing the calling script."""
    return Path(current_file).resolve().parent


def results_dir(current_file: str, subdir: str = "results/current") -> Path:
    """Create and return a results directory next to the scripts."""
    out = package_root(current_file) / subdir
    out.mkdir(parents=True, exist_ok=True)
    return out


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def environment_text() -> str:
    return "\n".join([
        f"created_utc = {datetime.now(timezone.utc).isoformat()}",
        f"python = {platform.python_version()}",
        f"platform = {platform.platform()}",
        f"processor = {platform.processor()}",
    ])
