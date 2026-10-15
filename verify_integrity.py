"""Verify every distributed file against SHA256SUMS, without third-party packages."""
from __future__ import annotations
import hashlib
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def verify() -> int:
    manifest=ROOT/'SHA256SUMS'
    entries=[]
    for line in manifest.read_text(encoding='utf-8').splitlines():
        if not line.strip():
            continue
        digest,name=line.split(maxsplit=1)
        path=(ROOT/name.lstrip('*')).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise AssertionError(f'Missing or invalid path: {name}')
        found=hashlib.sha256(path.read_bytes()).hexdigest()
        if found!=digest:
            raise AssertionError(f'SHA-256 mismatch: {name}')
        entries.append(path)
    if len(entries)!=len(set(entries)):
        raise AssertionError('Duplicate manifest entries.')
    if not entries:
        raise AssertionError('The integrity manifest is empty.')
    return len(entries)


if __name__=='__main__':
    print(f'PASS: {verify()} SHA-256 file checks.')
