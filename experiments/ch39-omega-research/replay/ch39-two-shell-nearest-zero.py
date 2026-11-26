"""Enclose the nearest other zero in the six normalized Wilson deficits on the two shells."""
from pathlib import Path
import runpy
from flint import acb, arb

base=runpy.run_path(str(Path(__file__).with_name("ch39-two-shell-finite-readout.py")))
q=base["q"]
book=base["book"]
def deficits(h):
    out=[]
    for (w,v),c in ((book["f"][:2],book["axis"]),(book["f"][2:],book["m"]["m"]["axis"])):
        a,b=(h*w).sin_cos()
        u,vv=(h*v).sin_cos()
        out.extend((a*a,u*u,1-(b*vv+c*a*u)**2))
    return out

def check():
    top=q(10)**8
    tol=q(".012")
    before,candidates=base["possibilities"](top,tol)
    target=deficits(base["h"])
    found={}
    for H in candidates:
        actual=deficits(H)
        if any(abs(v-a)>tol for v,a in zip(actual,target)):
            continue
        assert not acb(q(1)/2,H.lower()).zeta().contains(0)
        left=H.lower().zeta_nzeros().unique_fmpz()
        right=H.upper().zeta_nzeros().unique_fmpz()
        assert left is not None and right is not None
        for n in range(int(left)+1,int(right)+1):
            if n==719 or n in found:
                continue
            h=acb.zeta_zero(n).imag
            assert 0<h<top
            vals=deficits(h)
            diffs=[abs(v-a) for v,a in zip(vals,target)]
            d=max(x.lower() for x in diffs).union(max(x.upper() for x in diffs))
            found[n]=(d,h,diffs)
    assert found
    idx=min(found,key=lambda n:float(found[n][0].mid()))
    best,h,ds=found[idx]
    rest=min([tol.lower()]+[d.lower() for n,(d,_,_) in found.items() if n!=idx])
    assert best<rest and idx==136990954
    assert q(".00581637009466855")<best<q(".00581637009466857")
    print("candidate intervals",before,len(candidates),"other zeros in outer intervals",len(found))
    print("closest competing zero",idx,h.str(15),"maximum deficit deviation",best.str(15))
    print("signed-free deviations",*[v.str(9) for v in ds])
    print("lower bound on second competitor",rest.str(12))
if __name__=="__main__":
    check()
