"""Entire positive stiffness axis, CP ceiling, and Gaussian weights at zeta zero 68."""
from pathlib import Path
import runpy
from flint import arb

edge=runpy.run_path(str(Path(__file__).with_name("ch39-clock-edge-source.py")))
src=edge["src"]
cost,k,v,g=edge["evaluate"](edge["p"],edge["q"])
polys={f:src["loop_poly"](f,k,v) for f in src["all_faces"]()}
events=[((6,14,7),".001715",".001717"),((7,8,18),".11480",".11481"),
        ((8,12,7),".26467",".26469"),((7,8,18,12),".89412",".89414")]
cap=arb("3.75")
assert 1/g<cap


def off_cut(poly,a,b):
    z=src["value"](poly,a.union(b))
    if z.real>0 or not z.imag.contains(0):
        return
    assert b-a>arb("1e-11")
    mid=(a+b)/2
    off_cut(poly,a,mid)
    off_cut(poly,mid,b)


def isolate(poly,left,right):
    a,b=arb(left),arb(right)
    interval=a.union(b)
    derivative=[i*x for i,x in enumerate(poly) if i]
    assert src["value"](poly,interval).real<0
    assert not src["value"](derivative,interval).imag.contains(0)
    for _ in range(45):
        assert src["value"](poly,a).imag*src["value"](poly,b).imag<0
        mid=(a+b)/2
        if (src["value"](poly,mid).imag*src["value"](poly,a).imag)<0:
            b=mid
        else:
            a=mid
    return a.union(b)


def cp_ceiling(e,a,b,maximum):
    norms=[sum(abs(z)**2 for z in x) for x in (a,b)]
    t=arb(0).union(maximum)
    def divided_off(i,j):
        return (-e[i]*b[i]*a[j].conjugate()-a[i]*(e[j]*b[j]).conjugate()
                +t*norms[1]*a[i]*a[j].conjugate())
    numerator=maximum**3*abs((divided_off(0,1)*divided_off(1,2)*divided_off(2,0)).imag)
    correction=maximum*(norms[0]*norms[1]).sqrt()
    shift=2*max(abs(z).upper() for z in e)*correction+correction**2
    denominator=arb(1)
    for i in range(3):
        for j in range(i):
            gap=abs(abs(e[i])**2-abs(e[j])**2)-2*shift
            assert gap>0
            denominator*=gap
    return numerator/denominator


def probability(strength,count):
    reps=[(1,6),(0,7),(0,8),(1,8)]
    win=cost[reps[0]]
    return 1/sum((((1+strength*win)/(1+strength*cost[e])).log()*count).exp() for e in reps)


if __name__=="__main__":
    roots=[]
    for f,poly in polys.items():
        local=[(lo,hi) for a,lo,hi in events if a==f]
        start=arb(0)
        for lo,hi in local:
            off_cut(poly,start,arb(lo))
            roots.append(isolate(poly,lo,hi))
            start=arb(hi)
        off_cut(poly,start,cap)
    roots.sort(key=lambda z:float(z.mid()))
    examples=[0,".01",".18",".5","1"]
    observed=[src["degrees"](polys,arb(t)) for t in examples]
    assert observed==[(3,[3]),(2,[2]),(2,[1,2]),(2,[1]),(1,[1])]
    critical=roots[0]/(1-g*roots[0])
    j=cp_ceiling(*src["fourier_columns"](k,v),roots[0].upper())
    assert j<arb("3.2e-12")
    # For fixed S and positive integer multiplicity, limiting target probability is monotone in both.
    lo,hi=0,10000000
    while hi-lo>1:
        mid=(hi+lo)//2
        pr=probability(critical,mid)
        if pr>arb(".95"):
            hi=mid
        else:
            assert pr<arb(".95")
            lo=mid
    one=probability(critical,1)
    print("effective cuts",*[z.str(14) for z in roots])
    print("bare stiffness upper endpoint, J ceiling, one-copy orbit probability bound",
          critical.str(14),j.upper().str(12),one.upper().str(12))
    print("smallest copy count admitting 95 percent for some robust stiffness",hi)
