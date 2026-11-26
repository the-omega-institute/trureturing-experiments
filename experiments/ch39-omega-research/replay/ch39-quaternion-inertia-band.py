"""Enclose first 18 scalar levels for all longest-edge Dirichlet energies from 900 to infinity."""
from pathlib import Path
from fractions import Fraction
import runpy
from flint import arb,acb,acb_mat
root=Path(__file__).parent
m=runpy.run_path(str(root/'ch39-quaternion-scalar-levels.py'))
pencil=m['m']['m']['real_energy_matrix'].__globals__
ctxsrc=pencil['src']
ctxsrc['rho']=arb.pi()**2/(900*ctxsrc['theta'](804))
q=m['q']
def first18():
    result=[]
    for n in range(1,19):
        l,r=Fraction(1,8),Fraction(16)
        assert m['count'](q(l))==0 and m['count'](q(r))==18
        for _ in range(16):
            mid=(l+r)/2
            if m['count'](q(mid))<n:l=mid
            else:r=mid
        result.append((l,r))
    return result
lo=first18()
# Zero-inertia limit: the edge contributes its static spring matrix at every energy.
def static(s):
    mat=m['src']['matrix'](0,{})
    for j in range(19):mat[j,j]-=s
    return mat
m['count'].__globals__['m']['m']['real_energy_matrix']=static
hi=first18()
k=(16*arb.pi()**2/900).sqrt()
shift=4*arb('.01')*k/k.sin()
assert shift<arb('.0413')
boxes=[(q(a[0])-shift,q(b[1])+shift) for a,b in zip(lo,hi)]
assert boxes[0][0]>arb('1.6')
def fmt(v):
    return int(v.unique_fmpz())
if __name__=='__main__':
    path=root/'ch39-quaternion-band-boxes.txt'
    rows=[(fmt((x.lower()*10**6).floor()),fmt((y.upper()*10**6).ceil())) for x,y in boxes]
    path.write_text(''.join(f'{a} {b}\n' for a,b in rows))
    print(rows)
