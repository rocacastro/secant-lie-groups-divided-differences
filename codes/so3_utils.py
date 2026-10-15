"""Precision-safe SO(3) operations on the principal-logarithm domain."""
import mpmath as mp
try:
    from .lie_divdiff import norm
except ImportError:  # Support direct execution from the codes directory.
    from lie_divdiff import norm

def hat(v):
    x,y,z=v
    return mp.matrix([[0,-z,y],[z,0,-x],[-y,x,0]])
def vee(A):
    return mp.matrix([A[2,1],A[0,2],A[1,0]])
def exp(v):
    v=mp.matrix(v); t=norm(v); K=hat(v)
    if t==0:return mp.eye(3)
    return mp.eye(3)+(mp.sin(t)/t)*K+2*(mp.sin(t/2)/t)**2*(K*K)
def log(R):
    s=vee((R-R.T)/2); st=norm(s)
    ct=(mp.fsum(R[i,i] for i in range(3))-1)/2
    if st==0:
        if ct<=0: raise ValueError('Principal logarithm is not unique at angle pi.')
        return mp.zeros(3,1)
    return (mp.atan2(st,ct)/st)*s
