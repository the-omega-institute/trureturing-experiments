"""Propagate minimal positive Euclidean response realizations into the shell thresholds."""
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

directory=Path(__file__).parent
quad=run_path(str(directory/"ch39-kinetic-rh-gaussian-readout.py"))["gauss"]
base=run_path(str(directory/"ch39-finite-rh-propagating-references.py"))


def scalar_spectrum(mf,C,rho,p,modes):
    assert len(C)==len(rho) and all(u<v for u,v in zip(C,C[1:]))
    Z=(125/mf)**2
    spectrum=[(d,25) for d in C]
    for lam,n in modes:
        trace=lam+Z+sum(C)
        def secular(t):
            return t-lam*p-Z-lam*t*sum(r/(t-u) for r,u in zip(rho,C))
        # Each open interval between consecutive poles contains one simple
        # root; the last root is bounded above by the full positive trace.
        block=[]
        for left,right in zip([arb(0),*C],[*C,trace]):
            lo,hi=left.lower(),right.upper()
            for _ in range(72):
                mid=(lo+hi)/2
                assert all(not (mid-u).contains(0) for u in C)
                val=secular(mid)
                if val.contains(0):
                    break
                if val<0:
                    lo=mid
                else:
                    hi=mid
            root=lo.union(hi)
            assert left<root and root<right
            block.append(root)
            spectrum.append((root,n))
        assert (sum(block)-trace).contains(0)
    return spectrum


def main():
    ctx.prec,ctx.cap=3072,3
    a,b,c=[base["gram"](i,j) for i,j in ((1,0),(1,2),(0,2))]
    p=(1+2*a*b*c-a*a-b*b-c*c)/(1-c*c)
    modes=[(arb(6),7),(arb(4),6)]
    for k,n in ((5,3),(3,2)):
        modes+= [(3-arb(k).sqrt(),n),(3+arb(k).sqrt(),n)]
    d=[8*arb.pi()**2/arb(s)**2 for s in (".467453",".63811","1.0583")]
    K=5*d[0]-3*d[1]-2*d[2]
    F=arb(109)/15*(d[1]-d[2])-arb(23)/6*(d[0]-d[1])
    q=1000*(K/44).exp()
    for count,expected,low,high in [(2,71,"3.5515295e13","3.5515297e13"),(21,144,"7.4768884e13","7.4768887e13")]:
        z,weights,_=quad(count)
        av=sum(z)/count
        C=[av/v for v in z]
        rho=[(1-p)*w for w in weights]
        assert all(a<b for a,b in zip(C,C[1:])) and min(rho)>0
        for sign,mf in [(1,arb(low)),(-1,arb(high))]:
            y=q/mf
            assert y<arb(8).sqrt()
            spectrum=scalar_spectrum(mf,C,rho,p,modes)
            total=arb(0)
            active=0
            for v,n in spectrum:
                assert v>0 and mf*v.sqrt()>1000
                term=(y/v.sqrt()).log()
                if term>0:
                    total+=n*term
                    active+=n
                else:
                    assert term<0
            assert sum(n for _,n in spectrum)==23+48*count and active==expected
            assert sign*(22*total-15*F)>0
        print(count,"references per edge:",23+48*count,"extra scalar doublets,",expected,
              "open at Q =",q.str(14),"GeV; M_F bracket",low,high)


if __name__=="__main__":
    main()
