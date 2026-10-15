"""Verify reference data, exact bounds, metadata, and Bessel solver behavior.

This command checks archived SO(3)/SO(2) data without recomputing their long
trajectories. Pass --generated-dir after run_all.py to compare fresh results.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import re
import mpmath as mp
from codes.verify_semilocal import verify as verify_semilocal
from experiments.bessel.verify import verify as verify_bessel
from verify_integrity import verify as verify_integrity

ROOT=Path(__file__).resolve().parent
REFERENCE=ROOT/'results'/'reference'


def read_csv(path):
    with Path(path).open(newline='',encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def verify_metadata():
    version=(ROOT/'VERSION').read_text().strip()
    cff=(ROOT/'CITATION.cff').read_text(encoding='utf-8')
    cff_version=re.search(r'^version:\s*"?([^"\n]+)"?\s*$',cff,re.M).group(1)
    zenodo=json.loads((ROOT/'.zenodo.json').read_text(encoding='utf-8'))
    if version!=cff_version or version!=zenodo['version']:
        raise AssertionError('Version metadata are not synchronized.')
    if zenodo['license']!='MIT' or 'license: MIT' not in cff:
        raise AssertionError('License metadata are inconsistent.')
    orcids=re.findall(r'orcid: "https://orcid.org/([^"\n]+)"',cff)
    if orcids!=[x['orcid'] for x in zenodo['creators']]:
        raise AssertionError('Authorship metadata are inconsistent.')
    if len(orcids)!=3:
        raise AssertionError('Invalid author count.')
    title=re.search(r'^title:\s*"([^"\n]+)"\s*$',cff,re.M).group(1)
    if title!=zenodo['title']:
        raise AssertionError('Software titles are inconsistent.')
    return {'version':version,'author_count':len(orcids),'metadata_consistent':True}


def verify_so3():
    results={}
    with mp.workdps(100):
        for problem in ['alignment','averaging']:
            records=[]
            for dps in [1200,1400]:
                name=f'{problem}_red_{dps}'
                data=json.loads((REFERENCE/'so3'/f'{name}.json').read_text())
                rows=read_csv(REFERENCE/'so3'/f'{name}.csv')
                if len(rows)!=14 or [int(r['k']) for r in rows]!=list(range(14)):
                    raise AssertionError(f'Invalid trajectory indices: {name}')
                for row,record_row in zip(rows,data['rows']):
                    for key,value in row.items():
                        if str(record_row[key])!=value:
                            raise AssertionError(f'CSV/JSON disagreement: {name}, {key}')
                if data['updates']!=12 or data['function_calls']!=62 or data['terminal_index']!=13:
                    raise AssertionError(f'Incorrect update/evaluation counts: {name}')
                if mp.mpf(rows[-1]['residual_F'])>=mp.mpf('1e-700') or mp.mpf(rows[-2]['residual_F'])<=mp.mpf('1e-700'):
                    raise AssertionError(f'Incorrect stopping index: {name}')
                if mp.mpf(data['max_secant_defect'])>=mp.mpf('1e-1100'):
                    raise AssertionError(f'Large secant defect: {name}')
                if problem=='averaging':
                    # Directly check the objective entry and the limiting value.
                    if abs(mp.mpf(rows[1]['objective'])-mp.mpf('0.00546628663174382355357244854399974526819'))>mp.mpf('1e-40'):
                        raise AssertionError('Incorrect initial averaging objective.')
                    if abs(mp.mpf(rows[-1]['objective'])-mp.mpf(77)/30000)>mp.mpf('1e-50'):
                        raise AssertionError('Incorrect limiting averaging objective.')
                records.append(data)
            if any(a['residual_F']!=b['residual_F'] for a,b in zip(records[0]['rows'],records[1]['rows'])):
                raise AssertionError(f'Precision agreement failed: {problem}')
            results[problem]={'updates':12,'field_evaluations':62,'stored_residual_digits_match':50,
                              'final_residual':records[0]['rows'][-1]['residual_F']}
    return results


def compare_generated(output_dir):
    output=Path(output_dir)
    checked=[]
    for problem in ['alignment','averaging']:
        for dps in [1200,1400]:
            name=f'{problem}_red_{dps}'
            actual=json.loads((output/'so3'/f'{name}.json').read_text())
            expected=json.loads((REFERENCE/'so3'/f'{name}.json').read_text())
            for key in ['updates','function_calls','terminal_index','diagnostic_log_calls_for_objective']:
                if actual[key]!=expected[key]:
                    raise AssertionError(f'Generated count mismatch: {name}, {key}')
            if len(actual['rows'])!=len(expected['rows']):
                raise AssertionError(f'Generated length mismatch: {name}')
            for a,b in zip(actual['rows'],expected['rows']):
                for key in ['residual_F','rho_k','objective','angular_error','distance_to_X1']:
                    if a[key]!=b[key]:
                        raise AssertionError(f'Generated value mismatch: {name}, row {a["k"]}, {key}')
            if read_csv(output/'so3'/f'{name}.csv')!=read_csv(REFERENCE/'so3'/f'{name}.csv'):
                # Roundoff diagnostics (orthogonality and determinants) need not
                # agree bitwise; their scales are validated independently.
                with mp.workdps(100):
                    for row in actual['rows']:
                        if mp.mpf(row['orthogonality_F'])>=mp.mpf(10)**(-dps+10):
                            raise AssertionError(f'Large orthogonality defect: {name}')
                        if abs(mp.mpf(row['determinant'])-1)>mp.mpf('1e-45'):
                            raise AssertionError(f'Large determinant defect: {name}')
            checked.append(name)
    scalar=read_csv(output/'so2'/'so2_results.csv')
    if scalar!=read_csv(REFERENCE/'so2'/'so2_results.csv'):
        raise AssertionError('Scalar reference comparison failed.')
    checked.append('so2_results')
    return checked


def verify(generated_dir=None,run_solvers=True):
    result={'integrity_files_checked':verify_integrity(), 'metadata':verify_metadata(),
            'semilocal':verify_semilocal(), 'so3_reference':verify_so3(),
            'bessel':verify_bessel(run_solvers)}
    scalar=read_csv(REFERENCE/'so2'/'so2_results.csv')
    with mp.workdps(80):
        first=next(int(r['k']) for r in scalar if mp.mpf(r['residual_frobenius'])<=mp.mpf('1e-700'))
    if len(scalar)!=15 or first!=13:
        raise AssertionError('Scalar reference trajectory is inconsistent.')
    result['so2_reference']={'first_index_below_tolerance':first,'diagnostic_terminal_index':14}
    if generated_dir is not None:
        result['fresh_high_precision_comparisons']=compare_generated(generated_dir)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generated-dir',type=Path)
    parser.add_argument('--records-only',action='store_true',help='Skip Bessel solver execution.')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=verify(args.generated_dir,not args.records_only)
    if args.output:
        if args.output.resolve().is_relative_to(REFERENCE.resolve()):
            raise ValueError('Generated results must not overwrite reference data.')
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
