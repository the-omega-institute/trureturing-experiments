"""Global phase recovery from three spectral determinant readings, using python-flint."""
from fractions import Fraction
from pathlib import Path
import runpy
from flint import acb,acb_mat,arb,arb_mat
src=runpy.run_path(str(Path(__file__).with_name('ch39-arc-relative-determinant.py')))
def real_energy_matrix(s):
    m=acb_mat(19,19)
    for i,w in enumerate(src['V']):
        m[i,i]=w-s
    coeff={}
    for i,j,d,z in src['edges']:
        if d not in coeff:
            k=(s*src['rho']*src['theta'](d)).sqrt()
            sn,cs=k.sin_cos()
            assert not sn.contains(0)
            coeff[d]=k*cs/sn,k/sn
        diag,off=coeff[d]
        m[i,i]+=diag
        m[j,j]+=diag
        m[i,j]=-off*z
        m[j,i]=m[i,j].conjugate()
    return m
def coefficients(u,roots,modes):
    base=src['matrix'](u,{}) if u>=0 else real_energy_matrix(-u)
    clock_edges=[(i,j,d,base[i,j]) for i,j,d,_ in src['edges'] if d in (724,804)]
    norm=base.det().real
    assert not norm.contains(0)
    vals=[]
    for iz in range(7):
        for iw in range(7):
            M=acb_mat(base)
            for i,j,d,v in clock_edges:
                entry=v*roots[iz if d==724 else iw]
                M[i,j]=entry
                M[j,i]=entry.conjugate()
            val=M.det()
            assert val.imag.contains(0)
            vals.append(val.real/norm-1)
    return [2*sum(vals[7*i+j]*roots[(-a*i-b*j)%7].real
                  for i in range(7) for j in range(7))/49 for a,b in modes]
modes=[(a,b) for a in range(4) for b in range(-3,4)
       if (a>0 or b>0) and abs(a+b)<=3]
assert len(modes)==18
roots=[acb(0,2*arb.pi()*j/7).exp() for j in range(7)]

fund=[modes.index(k) for k in ((1,0),(0,1),(1,-1))]
extra=[j for j in range(18) if j not in fund]
# Integer monomials in x=cos(p), y=cos(q), z=cos(p-q).
def plus(p,q):
    r=p.copy()
    for k,val in q.items():r[k]=r.get(k,0)+val
    return {k:s for k,s in r.items() if s}
def scale(p,n):return {k:n*y for k,y in p.items()}
def mul(p,q):
    r={}
    for a,x in p.items():
        for b,y in q.items():
            k=tuple(i+j for i,j in zip(a,b))
            r[k]=r.get(k,0)+x*y
    return {k:s for k,s in r.items() if s}
one={(0,0,0):1};x={(1,0,0):1};y={(0,1,0):1};z={(0,0,1):1}
def cheb(v,n,second=False):
    p,q=one,scale(v,2 if second else 1)
    for _ in range(n):
        p,q=q,plus(scale(mul(v,q),2),scale(p,-1))
    return p
def cosine(a,b):
    C=mul(cheb(x,abs(a)),cheb(y,abs(b)))
    if a*b:
        U=mul(cheb(x,abs(a)-1,True),cheb(y,abs(b)-1,True))
        prod=plus(z,scale(mul(x,y),-1))
        C=plus(C,scale(mul(U,prod),-1 if a*b>0 else 1))
    return C
def lipschitz(mode):
    p=cosine(*mode)
    return sum(abs(val)*sum(k) for k,val in p.items())
L=[lipschitz(modes[j]) for j in extra]

def certify(shifts, contraction_limit, inverse_limit):
    mat=arb_mat([coefficients(arb(u),roots,modes) for u in shifts])
    base=arb_mat([[mat[i,j] for j in fund] for i in range(3)])
    remainder=arb_mat([[mat[i,j] for j in extra] for i in range(3)])
    inv=base.inv()
    mix=inv*remainder
    err=max(sum(abs(mix[i,j])*L[j] for j in range(len(extra))).upper() for i in range(3))
    amp=max(sum(abs(inv[i,j]) for j in range(3)).upper() for i in range(3))/(1-err)
    assert err<arb(contraction_limit) and amp<arb(inverse_limit)
    print(shifts,'contraction upper',err.str(12),'inverse sup norm upper',amp.upper().str(12))
    return mat
if __name__=='__main__':
    certify([15000,30000,100000],'.267','4.549e16')
    certify([-1000,-210000,-235000],'.176','3.161e11')
