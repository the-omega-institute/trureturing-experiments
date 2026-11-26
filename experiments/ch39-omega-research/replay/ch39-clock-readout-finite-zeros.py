"""Rule out every other zero below the given height by a sparse exact congruence and interval zeta."""
from fractions import Fraction
from pathlib import Path
import runpy
from flint import acb,arb,ctx

m=runpy.run_path(str(Path(__file__).with_name("ch39-clock-readout-resolution.py")))
model=m["m"]["m"]
sh=model["cur"]["sh"]
q,pi=m["q"],arb.pi()
w,v=[sh["theta"](d)-sh["theta"](244) for d in (724,804)]
ratio=v/w
height=acb.zeta_zero(719).imag
epsa=abs(m["hi"].sqrt().asin()-m["s"].sqrt().asin()).upper()
epsb=abs(m["lo"].sqrt().asin()-m["t"].sqrt().asin()).upper()
angular=(epsb+ratio*epsa)/pi

def congruence(top):
    K=int(((w*(top+height)+epsa)/pi).upper().ceil().unique_fmpz())
    # Only the proposal uses a rounded midpoint. Recompute the entire rational
    # approximation error against the enclosing frequency interval.
    base=10**32
    p=Fraction(int((base*ratio.mid()).floor().unique_fmpz()),base).limit_denominator(K)
    a,b=p.numerator,p.denominator
    radius=int((b*(angular+K*abs(ratio-q(a)/b))).upper().ceil().unique_fmpz())
    inv=pow(a,-1,b)
    cand=set()
    for residue in range(-radius,radius+1):
        k0=residue*inv % b
        for kk in range(k0-((K+k0)//b)*b,K+1,b):
            if abs(kk)>K:
                continue
            best=(ratio*kk).mid()
            ell=int((best+q(1)/2).floor().unique_fmpz())
            if not abs(ratio*kk-ell)>angular:
                cand.add((kk,ell))
    return cand

def windows(top):
    cand=congruence(top)
    intervals=[]
    for k,l in sorted(cand):
        if k==0:
            continue
        for sign in (-1,1):
            x=sign*height+(pi*k+arb(0,epsa))/w
            y=sign*height+(pi*l+arb(0,epsb))/v
            if x.overlaps(y):
                h=x.intersection(y)
                if not (h<0 or h>top):
                    intervals.append(h)
    return intervals

def exclude(h):
    pending=[(h,0)]
    leaves=0
    while pending:
        seg,depth=pending.pop()
        val=acb(q(1)/2,seg).zeta()
        if not val.contains(0):
            leaves+=1
            continue
        assert depth<20,seg
        mid=seg.mid()
        pending.extend((v,depth+1) for v in (seg.lower().union(mid),mid.union(seg.upper())))
    return leaves

def check(top):
    assert height-epsa/w>acb.zeta_zero(718).imag and height+epsa/w<acb.zeta_zero(720).imag
    intervals=windows(top)
    leaves=sum(exclude(h) for h in intervals)
    n=top.zeta_nzeros().unique_fmpz()
    assert n is not None
    print("angular halfwidths",epsa.str(10),epsb.str(10))
    print("maximum ordinate",top,"zero count",n,"other candidate boxes",len(intervals),"excluded pieces",leaves)

if __name__=="__main__":
    check(q(10)**8)
