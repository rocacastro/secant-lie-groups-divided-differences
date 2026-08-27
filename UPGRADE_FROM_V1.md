# Upgrading the GitHub repository from v1.0.0 to v2.0.0

The `SO(2)` scripts remain valid and should be retained. The two old `SO(3)`
scripts must be removed because they implement the component-wise divided
difference that is no longer used in the revised manuscript.

## Remove

```text
so3_secant_example2_restart_high_precision.py
so3_rotation_averaging_componentwise_alg3_high_precision.py
```

Also remove obsolete tracked files in `results/current/`; the revised repository
keeps generated results out of version control and stores manuscript reference
outputs in `results/reference/`.

## Add or replace

```text
README.md
CHANGELOG.md
CITATION.cff
.zenodo.json
.gitignore
requirements.txt
run_all.py
verify_reference.py
so3_rankone_divdiff_utils.py
so3_alignment_bd_rankone_high_precision.py
so3_rotation_averaging_r1_high_precision.py
results/README.md
results/reference/*
```

The three `SO(2)` scripts, `_repo_utils.py`, and `LICENSE` are retained.

After uploading the files, run:

```bash
python -m pip install -r requirements.txt
python run_all.py
python verify_reference.py
```

Then publish a GitHub release tagged `v2.0.0`. Zenodo should archive that release
and assign its DOI.
