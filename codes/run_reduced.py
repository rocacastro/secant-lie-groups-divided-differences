"""Compute the SO(3) alignment and rotation-averaging experiments.

    python run_reduced.py --problem alignment --dps 1200
    python run_reduced.py --problem averaging --dps 1200

The functions receive no derivative or structural information. Counts include
cached endpoint residuals and exclude the diagnostic objective J, whose cost is
recorded separately. No timing comparison is asserted.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import mpmath as mp
try:
    from .lie_divdiff import norm, reduced_divdiff
    from .so3_utils import exp, log
except ImportError:  # Support direct execution as a script.
    from lie_divdiff import norm, reduced_divdiff
    from so3_utils import exp, log


def run(problem='alignment',dps=1200,tol_digits=700,output_dir=None):
    if tol_digits < 1: raise ValueError('tol_digits must be positive.')
    if dps<tol_digits+150:raise ValueError('At least 150 guard digits are required.')
    mp.mp.dps=dps; calls=0; diag_calls=0
    if problem=='alignment':
        target=exp(mp.matrix([mp.pi/(4*mp.sqrt(3))]*3))
        b0=mp.matrix(list(map(mp.mpf,['-0.01309040431','-0.03172408338','0.00427275100'])))
        b1=mp.matrix(list(map(mp.mpf,['-0.03234918357','0.00015621005','-0.00459871328'])))
        rotations=None
    elif problem=='averaging':
        target=exp(mp.matrix(list(map(mp.mpf,['0.45','-0.25','0.35']))))
        rotations=[]
        for i,a in enumerate(['0.060','0.050','0.040']):
            for sign in [1,-1]:
                xi=mp.zeros(3,1);xi[i]=sign*mp.mpf(a);rotations.append(target*exp(xi))
        # The polar factor of the symmetric weighted average is exactly target.
        b0=mp.matrix(list(map(mp.mpf,['0.050','-0.040','0.030'])))
        b1=mp.matrix(list(map(mp.mpf,['-0.040','0.030','0.020'])))
    else:raise ValueError('Unknown problem.')
    def field(X):
        nonlocal calls
        calls+=1
        if problem=='alignment':
            value=log(X.T*target)
            if norm(value)>=mp.mpf(1)/5:
                raise RuntimeError('An alignment evaluation left the field domain.')
            return value
        out=mp.zeros(3,1)
        for R in rotations:out+=log(X.T*R)/6
        return out
    xs=[target*exp(b0),target*exp(b1)];fs=[field(X) for X in xs]
    res=[mp.sqrt(2)*norm(f) for f in fs];stats=[]
    tol=mp.mpf(10)**(-tol_digits)
    for k in range(1,40):
        before=calls
        L,v,Q,h=reduced_divdiff(field,xs[k-1],xs[k],fs[k-1],fs[k],exp,log)
        assert calls-before==4
        step=mp.lu_solve(L,-fs[k])
        if not all(mp.isfinite(x) for x in step) or norm(step)>=mp.pi:
            raise RuntimeError('A step is nonfinite or outside the logarithm branch.')
        stats.append({'k':k,'secant_defect':norm(L*v-(fs[k]-fs[k-1])),
                      'frame_defect':norm(Q.T*Q-mp.eye(3)),
                      'L_plus_I_frobenius':norm(L+mp.eye(3)),
                      'step_norm':norm(step),'h':h})
        Xnew=xs[k]*exp(step);xs.append(Xnew);fs.append(field(Xnew))
        res.append(mp.sqrt(2)*norm(fs[-1]))
        print(problem,dps,k+1,mp.nstr(res[-1],16),flush=True)
        if res[-1]<tol:break
    else:raise RuntimeError('Convergence tolerance not reached.')
    orders=[None]*len(xs)
    for k in range(1,len(xs)-1):
        den=mp.log(res[k]/res[k-1])
        if den:orders[k]=mp.log(res[k+1]/res[k])/den
    ns=lambda x: '--' if x is None else mp.nstr(x,50,min_fixed=0,max_fixed=0)
    rows=[]
    for k,X in enumerate(xs):
        err=norm(log(X.T*target));distance=norm(log(xs[1].T*X))
        if problem=='alignment' and k>=1:
            assert distance<mp.mpf('.06') and err<mp.mpf('.1')
        if problem=='averaging':
            J=mp.fsum(norm(log(X.T*R))**2 for R in rotations)/6;diag_calls+=6
        else:J=None
        rows.append({'k':k,'residual_F':ns(res[k]),'rho_k':ns(orders[k]),
                     'angular_error':ns(err),'distance_to_X1':ns(distance),
                     'objective':ns(J),'orthogonality_F':ns(norm(X.T*X-mp.eye(3))),
                     'determinant':ns(mp.det(X))})
    result={'problem':problem,'dps':dps,'tau':'0.1','tol_digits':tol_digits,
            'frame_rule':'signed Householder completion',
            'terminal_index':len(xs)-1,'updates':len(xs)-2,'function_calls':calls,
            'diagnostic_log_calls_for_objective':diag_calls,
            'mean_last_four_orders':ns(mp.fsum(orders[-5:-1])/4),
            'max_secant_defect':ns(max(t['secant_defect'] for t in stats)),
            'max_frame_defect':ns(max(t['frame_defect'] for t in stats)),
            'max_orthogonality_F':ns(max(norm(X.T*X-mp.eye(3)) for X in xs)),
            'max_L_plus_I_frobenius':ns(max(t['L_plus_I_frobenius'] for t in stats)),
            'first_step_norm':ns(stats[0]['step_norm']),
            'alpha':ns(norm(log(xs[0].T*xs[1]))),
            'rows':rows,'step_stats':[{key:(val if key=='k' else ns(val)) for key,val in s.items()} for s in stats]}
    assert calls==2+5*(len(xs)-2)
    output=(Path(output_dir) if output_dir is not None else
            Path(__file__).resolve().parents[1]/'results'/'generated'/'so3')
    reference = Path(__file__).resolve().parents[1]/'results'/'reference'
    if output.resolve().is_relative_to(reference.resolve()):
        raise ValueError('Generated results must not overwrite reference data.')
    output.mkdir(parents=True, exist_ok=True)
    name=f'{problem}_red_{dps}'
    (output/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    with (output/(name+'.csv')).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print('DONE',name,'calls',calls,'order',result['mean_last_four_orders'],flush=True)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--problem',choices=['alignment','averaging'],default='alignment');p.add_argument('--dps',type=int,default=1200);p.add_argument('--tol-digits',type=int,default=700)
    p.add_argument('--output-dir',type=Path,default=None)
    a=p.parse_args();run(a.problem,a.dps,a.tol_digits,a.output_dir)
