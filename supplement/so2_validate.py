"""High-precision scalar validation on SO(2).

The coordinate equation is g(theta)=2*sin(theta)+1/4=0. On the chosen local
angle branch, its Lie-group secant iteration is the scalar secant iteration.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import mpmath as mp


def run(output_dir=None, dps=1800) -> dict:
    if dps < 1800:
        raise ValueError('At least 1800 digits are required for the diagnostic iterate.')
    output = (Path(output_dir) if output_dir is not None else
              Path(__file__).resolve().parents[1]/'results'/'generated'/'so2')
    reference = Path(__file__).resolve().parents[1]/'results'/'reference'
    if output.resolve().is_relative_to(reference.resolve()):
        raise ValueError('Generated results must not overwrite reference data.')
    output.mkdir(parents=True, exist_ok=True)
    with mp.workdps(dps):
        def g(t):
            return 2*mp.sin(t)+mp.mpf(1)/4
        theta = [mp.mpf(0), mp.mpf(1)/7]
        # Index 14 is an additional diagnostic iterate used to compute rho_13.
        while len(theta) <= 14:
            a, b = theta[-2], theta[-1]
            ga, gb = g(a), g(b)
            if gb == ga:
                raise ArithmeticError('The scalar secant denominator vanished.')
            theta.append(b-gb*(b-a)/(gb-ga))
        residuals = [mp.sqrt(2)*abs(g(t)) for t in theta]
        rows = []
        for k, (t, r) in enumerate(zip(theta,residuals)):
            rho = ''
            if 1 <= k <= 13:
                rho = mp.nstr(mp.log(residuals[k+1]/r)/mp.log(r/residuals[k-1]),16)
            rows.append({'k':k,'theta':mp.nstr(t,40),
                         'residual_frobenius':mp.nstr(r,16),'rho_k':rho})
        first = next(k for k,r in enumerate(residuals) if r <= mp.mpf('1e-700'))
        orders = [mp.log(residuals[k+1]/residuals[k])/mp.log(residuals[k]/residuals[k-1])
                  for k in range(9,14)]
        result = {'dps':dps,'tolerance':'1e-700','first_index_below_tolerance':first,
                  'updates_to_tolerance':first-1,'diagnostic_terminal_index':14,
                  'mean_orders_9_to_13':mp.nstr(mp.fsum(orders)/len(orders),16),
                  'known_root':mp.nstr(-mp.asin(mp.mpf(1)/8),40),
                  'note':'Index 14 is used only to compute the last order indicator.'}
    with (output/'so2_results.csv').open('w',newline='',encoding='utf-8') as stream:
        writer = csv.DictWriter(stream,fieldnames=rows[0].keys())
        writer.writeheader(); writer.writerows(rows)
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--dps',type=int,default=1800)
    args=parser.parse_args()
    run(args.output_dir,args.dps)
