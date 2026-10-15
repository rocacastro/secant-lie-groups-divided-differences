"""Bessel field, cached secant and Newton solvers, and paired timing experiment.

Default outputs are written under results/generated/bessel; reference records
are never overwritten. Run with --help for the output-directory option.
"""
from __future__ import annotations
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('MKL_NUM_THREADS','1')
os.environ.setdefault('NUMEXPR_NUM_THREADS','1')

from pathlib import Path
from time import perf_counter_ns
import json, hashlib, platform, csv
import numpy as np
import scipy
from scipy.linalg import expm, logm
from scipy.special import jv
from threadpoolctl import threadpool_limits, threadpool_info

ROOT=Path(__file__).resolve().parents[2] / 'results' / 'generated' / 'bessel'
BASE=Path(__file__).resolve().parent / 'problem_dense_original.npz'
P=100
ALPHA=0.15
SCALE=0.4
OFFSET=3.0
RADIUS=0.10
DELTA=0.01
TOL=1e-12
MAX_UPDATES=20
BLOCKS=51
REPS=5
WARMUP=20

class Field:
    def __init__(self, data):
        self.C=np.ascontiguousarray(data['C']); self.D=np.ascontiguousarray(data['D']); self.A=np.ascontiguousarray(data['A'])
        self.Rstar=np.ascontiguousarray(data['Rstar']); Xorig=np.ascontiguousarray(data['Xstart'])
        self.n=self.A.shape[0]; self.i,self.j=np.triu_indices(self.n,1); self.m=len(self.i); self.Avec=self.A[self.i,self.j].copy()
        k=np.arange(P,dtype=float)
        self.nu=k.copy()
        self.kappa=(k+1.0)/3.0
        self.weight=(1.0+0.25*np.cos(0.7*k))/P
        ir,jr=self.i[:,None],self.j[:,None]; ac,bc=self.i[None,:],self.j[None,:]
        self.e1=(jr==bc); self.e2=(jr==ac); self.e3=(ir==ac); self.e4=(ir==bc)
        # Deterministic start: use the supplied reference direction at fixed radius.
        V=logm(self.Rstar.T@Xorig).real
        u=V[self.i,self.j]; u=u/np.linalg.norm(u)
        self.Xstart=self.Rstar@expm(self.hat(RADIUS*u))
        gR=self.g(self.Rstar); assert np.linalg.norm(gR)==0.0
        wR=OFFSET+SCALE*(self.D@gR)
        phiR,_=self.phi_vals(wR)
        self.K=self.C@gR+ALPHA*phiR
        self.reset()
    def reset(self):
        self.F_calls=self.J_calls=self.sec_builds=self.exp_calls=0
        self.bessel_F_batches=self.bessel_J_batches=0; self.fallbacks=0; self.secant_defects=[]
    def hat(self,v):
        V=np.zeros((self.n,self.n),dtype=float); V[self.i,self.j]=v; V[self.j,self.i]=-v; return V
    def g(self,X):
        return X[self.i,self.j]-X[self.j,self.i]+self.Avec
    def G(self,Z):
        ir,jr=self.i[:,None],self.j[:,None]; ac,bc=self.i[None,:],self.j[None,:]
        return (Z[ir,ac]*self.e1-Z[ir,bc]*self.e2+Z[jr,bc]*self.e3-Z[jr,ac]*self.e4)
    def phi_vals(self,w):
        x=self.kappa[:,None]*w[None,:]
        vals=jv(self.nu[:,None],x)
        return self.weight@vals, vals
    def eval(self,X):
        self.F_calls+=1; self.bessel_F_batches+=1
        g=self.g(X); w=OFFSET+SCALE*(self.D@g); phi,vals=self.phi_vals(w)
        f=self.C@g+ALPHA*phi-self.K
        return f,w,phi,vals
    def dphi_from_cache(self,w,vals):
        # Optimized exact derivative, reusing J_nu from F:
        # J'_nu(z)=nu/z J_nu(z)-J_{nu+1}(z).  Only one new Bessel batch.
        x=self.kappa[:,None]*w[None,:]
        self.bessel_J_batches+=1
        jp=(self.nu[:,None]/x)*vals-jv((self.nu+1.0)[:,None],x)
        return self.weight@(self.kappa[:,None]*jp)
    def jacobian(self,X,w,vals):
        self.J_calls+=1
        dp=self.dphi_from_cache(w,vals)
        S=self.C+(ALPHA*SCALE*dp)[:,None]*self.D
        return S@self.G(X)
    def phi_derivative_direct(self,w):
        x=self.kappa[:,None]*w[None,:]
        vals=jv(self.nu[:,None],x)
        jp=(self.nu[:,None]/x)*vals-jv((self.nu+1.0)[:,None],x)
        return self.weight@(self.kappa[:,None]*jp)
    def alg_secant(self,X0,X1,f0,w0,p0,f1,w1,p1,v):
        self.sec_builds+=1
        r=float(np.linalg.norm(v)); q=v/r; dw=w1-w0
        thresh=1e-12*(1.0+np.maximum(np.abs(w0),np.abs(w1)))
        slope=np.empty_like(dw); mask=np.abs(dw)>thresh
        slope[mask]=(p1[mask]-p0[mask])/dw[mask]
        if np.any(~mask):
            self.fallbacks+=int(np.sum(~mask))
            slope[~mask]=self.phi_derivative_direct(0.5*(w0[~mask]+w1[~mask]))
        S=self.C+(ALPHA*SCALE*slope)[:,None]*self.D
        B=S@self.G(0.5*(X0+X1))
        L=B+((f1-f0)/r-B@q)[:,None]*q[None,:]
        self.secant_defects.append(float(np.linalg.norm(L@v-(f1-f0))/(1+np.linalg.norm(f1-f0))))
        return L
    def step(self,X,s):
        self.exp_calls+=1
        return X@expm(self.hat(s))

def solve_newton(p):
    p.reset(); X=p.Xstart.copy(); f,w,ph,vals=p.eval(X); tr=[float(np.linalg.norm(f))]
    for _ in range(MAX_UPDATES):
        if tr[-1]<=TOL: break
        J=p.jacobian(X,w,vals); s=np.linalg.solve(J,-f)
        if not np.all(np.isfinite(s)) or np.linalg.norm(s)>=np.pi: raise RuntimeError('Newton branch failure')
        X=p.step(X,s); f,w,ph,vals=p.eval(X); tr.append(float(np.linalg.norm(f)))
    if tr[-1]>TOL: raise RuntimeError('Newton failed')
    return X,tr,{k:getattr(p,k) for k in ['F_calls','J_calls','sec_builds','exp_calls','bessel_F_batches','bessel_J_batches','fallbacks']}

def solve_secant(p):
    p.reset(); v=DELTA*np.ones(p.m)/np.sqrt(p.m)
    X0=p.step(p.Xstart,-v); f0,w0,p0,_=p.eval(X0)
    X=p.Xstart.copy(); f,w,ph,_=p.eval(X); tr=[float(np.linalg.norm(f))]
    for _ in range(MAX_UPDATES):
        if tr[-1]<=TOL: break
        L=p.alg_secant(X0,X,f0,w0,p0,f,w,ph,v); s=np.linalg.solve(L,-f)
        if not np.all(np.isfinite(s)) or np.linalg.norm(s)>=np.pi: raise RuntimeError('Secant branch failure')
        Y=p.step(X,s); fy,wy,py,_=p.eval(Y)
        X0,f0,w0,p0,v=X,f,w,ph,s
        X,f,w,ph=Y,fy,wy,py; tr.append(float(np.linalg.norm(f)))
    if tr[-1]>TOL: raise RuntimeError('Secant failed')
    return X,tr,{k:getattr(p,k) for k in ['F_calls','J_calls','sec_builds','exp_calls','bessel_F_batches','bessel_J_batches','fallbacks']},max(p.secant_defects or [0.0])

def validate_jacobian(p):
    X=p.Xstart; f,w,ph,vals=p.eval(X); J=p.jacobian(X,w,vals)
    rng=np.random.default_rng(20261004); errs=[]
    for _ in range(5):
        u=rng.normal(size=p.m); u/=np.linalg.norm(u); h=1e-6
        fp,*_=p.eval(X@expm(p.hat(h*u))); fm,*_=p.eval(X@expm(p.hat(-h*u)))
        fd=(fp-fm)/(2*h); errs.append(np.linalg.norm(J@u-fd)/(1+np.linalg.norm(fd)))
    return max(errs)

def mode_contributions(p):
    out={}
    for label,X in [('root',p.Rstar),('start',p.Xstart)]:
        g=p.g(X); w=OFFSET+SCALE*(p.D@g); x=p.kappa[:,None]*w[None,:]
        c=np.abs(p.weight[:,None]*jv(p.nu[:,None],x)); rowmax=c.max(axis=1)
        out[label]={'w_min':float(w.min()),'w_max':float(w.max()),'min_row_max':float(rowmax.min()),'median_row_max':float(np.median(rowmax)),'max_row_max':float(rowmax.max()),'rows_gt_1e-5':int(np.sum(rowmax>1e-5))}
    return out

def time_one(fn,p,reps):
    t0=perf_counter_ns()
    for _ in range(reps): fn(p)
    return (perf_counter_ns()-t0)/reps/1e6

def main():
    global ROOT
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT)
    args = parser.parse_args()
    ROOT = args.output_dir
    reference = Path(__file__).resolve().parents[2]/'results'/'reference'
    if ROOT.resolve().is_relative_to(reference.resolve()):
        raise ValueError('Generated results must not overwrite reference data.')
    ROOT.mkdir(parents=True, exist_ok=True)
    with np.load(BASE, allow_pickle=False) as data:
        p=Field(data)
    with threadpool_limits(limits=1):
        Xn,tn,cn=solve_newton(p); Xs,ts,cs,defect=solve_secant(p)
        jac_err=validate_jacobian(p); modes=mode_contributions(p)
        assert cn['F_calls']==len(tn) and cn['J_calls']==len(tn)-1
        assert cn['bessel_F_batches']==cn['F_calls'] and cn['bessel_J_batches']==cn['J_calls']
        assert cs['F_calls']==len(ts)+1 and cs['J_calls']==0 and cs['bessel_J_batches']==0 and cs['fallbacks']==0
        for _ in range(WARMUP): solve_newton(p); solve_secant(p)
        rows=[]
        for b in range(BLOCKS):
            order=[('newton',solve_newton),('secant',solve_secant)] if b%2==0 else [('secant',solve_secant),('newton',solve_newton)]
            row={'block':b}
            for name,fn in order: row[name]=time_one(fn,p,REPS)
            row['ratio']=row['secant']/row['newton']; rows.append(row)
    def st(key):
        a=np.array([r[key] for r in rows]); return {'median_ms':float(np.median(a)),'q1_ms':float(np.quantile(a,.25)),'q3_ms':float(np.quantile(a,.75)),'min_ms':float(a.min()),'max_ms':float(a.max())}
    ratio=np.array([r['ratio'] for r in rows])
    out={
      'definition':{'n':p.n,'m':p.m,'P':P,'alpha':ALPHA,'scale':SCALE,'offset':OFFSET,'radius':RADIUS,'delta':DELTA,'tol':TOL,
                    'nu':'0,...,P-1','kappa':'(nu+1)/3','weights':'(1+0.25*cos(0.7*nu))/P'},
      'protocol':{'blocks':BLOCKS,'reps_per_block':REPS,'warmup_pairs':WARMUP,'blas_threads':1,'base_sha256':hashlib.sha256(BASE.read_bytes()).hexdigest()},
      'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform(),'threadpools':threadpool_info()},
      'validation':{
        'newton':{'updates':len(tn)-1,'residuals':tn,'root_distance_fro':float(np.linalg.norm(Xn-p.Rstar)),**cn},
        'secant':{'updates':len(ts)-1,'residuals':ts,'root_distance_fro':float(np.linalg.norm(Xs-p.Rstar)),'max_secant_defect':defect,**cs},
        'final_solution_distance':float(np.linalg.norm(Xn-Xs)),'max_directional_jacobian_rel_error':jac_err,'mode_contributions':modes},
      'timing':{'newton':st('newton'),'secant':st('secant'),'median_ratio_secant_over_newton':float(np.median(ratio)),'q1_ratio':float(np.quantile(ratio,.25)),'q3_ratio':float(np.quantile(ratio,.75)),'secant_faster_blocks':int(np.sum(ratio<1)),'blocks':BLOCKS},
      'blocks':rows}
    (ROOT/'Bessel_field_benchmark.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    with (ROOT/'Bessel_field_blocks.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['block','newton','secant','ratio']);w.writeheader();w.writerows(rows)
    print(json.dumps({k:v for k,v in out.items() if k!='blocks'},indent=2))

if __name__=='__main__': main()
