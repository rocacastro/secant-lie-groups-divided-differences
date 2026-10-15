"""Check the rational SO(3) semilocal bounds using exact arithmetic.

The uniform analytic estimate and its domain assumptions are given in
`docs/semilocal_bounds.md`. This script checks their rational consequences,
not the validity of an inequality by sampling points.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

B0 = tuple(Q(s) for s in ('-0.01309040431', '-0.03172408338', '0.00427275100'))
B1 = tuple(Q(s) for s in ('-0.03234918357', '0.00015621005', '-0.00459871328'))


def verify() -> dict:
    radius = Q(1, 10)
    probe_radius = Q(3, 25)
    delta = lambda t: t/2 + t*t/10
    d0, d1 = delta(radius), delta(probe_radius)
    sigma_sq = d0*d0 + 2*d1*d1
    sigma = Q(3, 25)
    omega = 2*sigma
    beta = 1/(1-sigma)
    eta = Q(3, 80)
    q = beta*omega
    a = b = c = q/(1-q)
    R = (a*b/(1-c) + a + 1)*eta
    norm0_sq = sum(t*t for t in B0)
    norm1_sq = sum(t*t for t in B1)
    cross = (B0[1]*B1[2]-B0[2]*B1[1],
             B0[2]*B1[0]-B0[0]*B1[2],
             B0[0]*B1[1]-B0[1]*B1[0])
    averaging_bound = sum(delta(t) for t in (Q(6,100), Q(5,100), Q(4,100)))/3
    checks = {
        'gamma_bound': 1/(12*(1-probe_radius**2/24)) < Q(1,10),
        'uniform_operator_bound': sigma_sq < sigma*sigma,
        'initial_X0_in_domain': norm0_sq < Q(35,1000)**2 < radius**2,
        'initial_X1_in_domain': norm1_sq < Q(33,1000)**2 < radius**2,
        'inverse_bound': sigma < 1,
        'first_step_bound': beta*Q(33,1000) == eta,
        'constant_omega': omega == Q(6,25),
        'unique_radius_root': R == Q(3,50),
        'beta_times_omega': q == Q(3,11) and q < 1,
        'contraction_factor': c == Q(3,8) and c < 1,
        'closed_ball_in_domain': Q(33,1000)+R < radius,
        # pi > 3, so these rational inequalities place the increments in W.
        'increment_radius_below_pi': R < 3,
        'pairwise_angle_below_pi': 2*radius < 3,
        'nonparallel_logarithmic_errors': any(t != 0 for t in cross),
        'averaging_derivative_invertible': averaging_bound == Q(7577,300000) < 1,
        'averaging_objective_at_root': (Q(6,100)**2+Q(5,100)**2+Q(4,100)**2)/3 == Q(77,30000),
    }
    if not all(checks.values()):
        raise ArithmeticError(f'Failed rational checks: {checks}')
    exact = dict(domain_radius=radius, probe_radius=probe_radius,
                 delta_segment=d0, delta_probe=d1, sigma_squared_bound=sigma_sq,
                 sigma=sigma, omega=omega, beta=beta, eta=eta, beta_omega=q,
                 a=a, b=b, c=c, R=R, norm_b0_squared=norm0_sq, norm_b1_squared=norm1_sq,
                 containment_margin_lower=radius-Q(33,1000)-R,
                 averaging_derivative_bound=averaging_bound,
                 averaging_objective_at_root=Q(77,30000))
    return {'checks':checks, 'constants_exact':{k:str(v) for k,v in exact.items()},
            'initial_b0':[str(t) for t in B0], 'initial_b1':[str(t) for t in B1],
            'scope':'Exact rational inequalities under the analytic bounds in docs/semilocal_bounds.md.'}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify()
    if args.output:
        reference = Path(__file__).resolve().parents[1]/'results'/'reference'
        if args.output.resolve().is_relative_to(reference.resolve()):
            raise ValueError('Generated results must not overwrite reference data.')
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    for name in result['checks']:
        print(f'PASS {name}')
    print('beta=25/22; eta=3/80; omega=6/25; a=b=c=3/8; R=3/50.')


if __name__ == '__main__':
    main()
