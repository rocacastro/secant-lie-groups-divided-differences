"""Verify the Bessel reference records and optionally run both solvers.

Statistics are recalculated from the 51 recorded blocks. Solver checks do not
rerun the timing experiment or replace its reference values.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path
import numpy as np

REPOSITORY = Path(__file__).resolve().parents[2]
REFERENCE = REPOSITORY/'results'/'reference'/'bessel'


def verify(run_solvers: bool = True) -> dict:
    record = json.loads((REFERENCE/'benchmark.json').read_text(encoding='utf-8'))
    hardware = json.loads((REFERENCE/'hardware.json').read_text(encoding='utf-8'))
    with (REFERENCE/'blocks.csv').open(newline='',encoding='utf-8') as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != record['protocol']['blocks'] or len(rows) != 51:
        raise AssertionError('The reference must contain 51 paired blocks.')
    if [int(r['block']) for r in rows] != list(range(51)):
        raise AssertionError('Reference block indices are not consecutive.')
    newton = np.array([float(r['newton']) for r in rows])
    secant = np.array([float(r['secant']) for r in rows])
    if not np.all(np.isfinite(newton)) or not np.all(newton > 0) or not np.all(np.isfinite(secant)) or not np.all(secant > 0):
        raise AssertionError('Timing values must be finite and positive.')
    ratios=secant/newton
    np.testing.assert_allclose(ratios,[float(r['ratio']) for r in rows],rtol=1e-14,atol=0)
    for i,row in enumerate(record['blocks']):
        if row['block']!=i or row['newton']!=newton[i] or row['secant']!=secant[i]:
            raise AssertionError('CSV and JSON timing records disagree.')
    for name,values in [('newton',newton),('secant',secant)]:
        stats=record['timing'][name]
        np.testing.assert_allclose(np.quantile(values,[.25,.5,.75]),
            [stats['q1_ms'],stats['median_ms'],stats['q3_ms']],rtol=1e-14,atol=0)
    paired=np.quantile(ratios,[.25,.5,.75])
    np.testing.assert_allclose(paired,[record['timing']['q1_ratio'],
        record['timing']['median_ratio_secant_over_newton'],record['timing']['q3_ratio']],rtol=1e-14,atol=0)
    if int((ratios<1).sum())!=record['timing']['secant_faster_blocks'] or int((ratios<1).sum())!=50:
        raise AssertionError('The paired timing count is inconsistent.')
    problem_path=REPOSITORY/'experiments'/'bessel'/'problem_dense_original.npz'
    digest=hashlib.sha256(problem_path.read_bytes()).hexdigest()
    if digest != record['protocol']['base_sha256']:
        raise AssertionError('The Bessel input does not match the recorded experiment.')
    result={'reference_statistics_verified':True,'input_sha256_verified':True,
        'hardware_record_present':bool(hardware),'median_paired_ratio':float(paired[1]),
        'paired_ratio_quartiles':[float(paired[0]),float(paired[2])],
        'paired_reduction_percent':float(100*(1-paired[1])),
        'secant_faster_blocks':int((ratios<1).sum()),'timing_rerun_performed':False,
        'solvers_executed':run_solvers}
    if run_solvers:
        import scipy
        from threadpoolctl import threadpool_limits,threadpool_info
        try:
            from . import benchmark as kernel
        except ImportError:
            import benchmark as kernel
        raw_jv=kernel.jv
        jv_calls=[0]
        def counted_jv(*args,**kwargs):
            jv_calls[0]+=1
            return raw_jv(*args,**kwargs)
        assembly_checks=[]
        kernel.jv=counted_jv
        try:
            with threadpool_limits(limits=1):
                with np.load(problem_path,allow_pickle=False) as data:
                    field=kernel.Field(data)
                xn,tn,cn=kernel.solve_newton(field)
                newton_jv_calls=jv_calls[0]
                # Count solve calls only; setup calls are not charged to either method.
                setup_calls=1
                if newton_jv_calls-setup_calls != 7:
                    raise AssertionError('Unexpected number of Bessel calls in Newton.')
                base_assembly=field.alg_secant
                keys=['F_calls','J_calls','bessel_F_batches','bessel_J_batches','fallbacks']
                def checked_assembly(*args,**kwargs):
                    before=tuple(getattr(field,k) for k in keys)+(jv_calls[0],)
                    answer=base_assembly(*args,**kwargs)
                    after=tuple(getattr(field,k) for k in keys)+(jv_calls[0],)
                    assembly_checks.append(before==after)
                    if before!=after:
                        raise AssertionError('Unexpected field, Jacobian, or Bessel evaluation in an assembly.')
                    return answer
                field.alg_secant=checked_assembly
                xs,ts,cs,defect=kernel.solve_secant(field)
                secant_jv_calls=jv_calls[0]-newton_jv_calls
                if secant_jv_calls != 6:
                    raise AssertionError('Unexpected number of Bessel calls in the secant solve.')
                jac_error=kernel.validate_jacobian(field)
                pools=threadpool_info()
        finally:
            kernel.jv=raw_jv
        if len(tn)-1!=3 or len(ts)-1!=4 or len(assembly_checks)!=4 or not all(assembly_checks):
            raise AssertionError('Solver update or assembly counts differ from the reference.')
        for method,counts in [('newton',cn),('secant',cs)]:
            for key,value in counts.items():
                if value!=record['validation'][method][key]:
                    raise AssertionError((method,key,value))
        if not (tn[-1]<=kernel.TOL and ts[-1]<=kernel.TOL and defect<1e-13 and jac_error<1e-7):
            raise AssertionError('A numerical accuracy check failed.')
        if np.linalg.norm(xn-field.Rstar)>=1e-10 or np.linalg.norm(xs-field.Rstar)>=1e-10:
            raise AssertionError('A solver did not reach the known root.')
        # Residuals at roundoff are checked by tolerance, not bitwise comparison.
        result.update({'counts_reproduced':True,'all_four_assemblies_without_extra_calls':True,
            'newton_check':{'residuals':tn,'counts':cn,'bessel_calls_in_solve':7},
            'secant_check':{'residuals':ts,'counts':cs,'bessel_calls_in_solve':6,'secant_defect':float(defect)},
            'jacobian_directional_discrepancy':float(jac_error),
            'validation_environment':{'python':platform.python_version(),'numpy':np.__version__,
                'scipy':scipy.__version__,'platform':platform.platform(),'threadpools':pools}})
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stats-only',action='store_true')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=verify(not args.stats_only)
    if args.output:
        if args.output.resolve().is_relative_to((REPOSITORY/'results'/'reference').resolve()):
            raise ValueError('Generated results must not overwrite reference data.')
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
