"""Run reproducibility checks and high-precision experiments.

The default recomputes SO(3) at 1200 and 1400 digits, SO(2) at 1800 digits,
and both Bessel solvers without timing. --quick skips high-precision reruns.
Use --benchmark only when a fresh timing experiment is intended.
"""
from __future__ import annotations
import argparse
import json
import platform
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent


def execute(arguments):
    print('+ '+ ' '.join(str(x) for x in [sys.executable,*arguments]),flush=True)
    subprocess.run([sys.executable,*map(str,arguments)],cwd=ROOT,check=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--quick',action='store_true')
    parser.add_argument('--benchmark',action='store_true')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'results'/'generated')
    args=parser.parse_args()
    output=args.output_dir.resolve()
    if output.is_relative_to((ROOT/'results'/'reference').resolve()):
        parser.error('Generated results must not overwrite reference data.')
    output.mkdir(parents=True,exist_ok=True)
    execute(['verify_integrity.py'])
    execute(['-m','unittest','discover','-s','tests','-v'])
    execute(['codes/verify_semilocal.py','--output',output/'semilocal.json'])
    if not args.quick:
        for problem in ['alignment','averaging']:
            for dps in [1200,1400]:
                execute(['-m','codes.run_reduced','--problem',problem,'--dps',dps,
                         '--output-dir',output/'so3'])
        execute(['-m','supplement.so2_validate','--output-dir',output/'so2'])
    verify_args=['verify_reference.py','--output',output/'verification.json']
    if not args.quick:
        verify_args+=['--generated-dir',output]
    execute(verify_args)
    if args.benchmark:
        execute(['-m','experiments.bessel.benchmark','--output-dir',output/'bessel'])
    import mpmath,numpy,scipy,threadpoolctl
    environment={'python':platform.python_version(),'platform':platform.platform(),
                 'mpmath':mpmath.__version__,'numpy':numpy.__version__,
                 'scipy':scipy.__version__,'threadpoolctl':threadpoolctl.__version__,
                 'high_precision_experiments_recomputed':not args.quick,
                 'timing_experiment_recomputed':args.benchmark}
    (output/'environment.json').write_text(json.dumps(environment,indent=2)+'\n',encoding='utf-8')
    print(f'All requested checks passed. Outputs: {output}',flush=True)


if __name__=='__main__':
    main()
