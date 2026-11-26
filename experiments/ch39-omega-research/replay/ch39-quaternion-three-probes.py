"""Three fixed spectral determinants distinguish simultaneous conjugacy on the full SU(2) pair space."""
import runpy
from pathlib import Path
from flint import arb,arb_mat
root=Path(__file__).parent
m=runpy.run_path(str(root/'ch39-quaternion-spectral-polynomial.py'))
one=(0,0,0)
axes=[(1,0,0),(0,1,0),(0,0,1)]
def certify(shifts):
    polys=[m['polynomial'](u) for u in shifts]
    keys=sorted(set(polys[0])-{one,*axes})
    A=arb_mat([[p[i] for i in axes] for p in polys])
    B=arb_mat([[p[i] for i in keys] for p in polys])
    inv=A.inv()
    rem=inv*B
    L=max(sum(abs(rem[i,j])*sum(k) for j,k in enumerate(keys)).upper() for i in range(3))
    assert L<1
    amp=max(sum(abs(inv[i,j]) for j in range(3)).upper() for i in range(3))/(1-L)
    return L,amp
if __name__=='__main__':
    for shifts, target, amp_limit in (([15000,20000,60000],'.400','1.187e16'),([-1000,-204000,-215000],'.331','5.196e10')):
        L,amp=certify(shifts)
        assert L<arb(target) and amp<arb(amp_limit)
        print(shifts,L.str(12),amp.str(12))
