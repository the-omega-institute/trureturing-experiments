"""Last common denominator for which fully relaxed finite reference angles change the winning zero."""
from itertools import product
from functools import lru_cache
from pathlib import Path
from runpy import run_path
from flint import arb,acb

home=Path(__file__).resolve().parent
model=run_path(str(home/"ch39-shell-continuous-penalty-switch.py"))
m=model["m"]
pi=arb.pi()
lam=m["penalty"]
relaxed={25:model["e25"]+lam*model["h25"]**2,255:model["e255"]+lam*model["h255"]**2}
obs={n:m["state"](acb.zeta_zero(n).imag) for n in relaxed}
heights={n:acb.zeta_zero(n).imag for n in relaxed}
center={25:(".452464611264",".562344350051",".117109250530"),
        255:(".162097469513",".546059802744",".132682748073")}


def atan_ref(n):
    h=heights[n]
    return [arb(v) for v in center[n]]+[h*f/pi for f in m["m"]["shells"][867][2]]


def tables(n,order):
    @lru_cache(None)
    def trig(i):
        return (pi*(i%order)/order).sin_cos()
    def residual(r,i,j):
        c=m["m"]["shells"][r][1]
        s,x=trig(i); t,y=trig(j)
        v=obs[n][r]
        return (s*s-v[0])**2+(t*t-v[1])**2+(1-(x*y+c*s*t)**2-v[2])**2
    return trig,residual


def trial(n,order,neighbors):
    _,loss=tables(n,order)
    indices=[round(float(z.mid())*order) for z in atan_ref(n)]
    values=[tuple(range(j-1,j+2)) if neighbors else (j,) for j in indices]
    small=min((loss(17,p,q)+loss(26,a,q)+loss(29,a,p)).upper() for p,q,a in product(*values[:3]))
    big=min(loss(867,i,j).upper() for i,j in product(*values[3:]))
    return small+big+lam*heights[n]**2


def lower(n,order,upper):
    trig,loss=tables(n,order)
    ys=obs[n]
    hcost=lam*heights[n]**2
    # Discard a coordinate only when its contribution alone precludes an improvement.
    bound=m["paired_lower"](ys)
    means=[(ys[17][0]+ys[29][1])/2,(ys[17][1]+ys[26][1])/2,(ys[26][0]+ys[29][0])/2]
    def options(v,two):
        return [j for j in range(order) if not (2 if two else 1)*(trig(j)[0]**2-v)**2+bound+hcost>upper]
    smallopts=[options(v,True) for v in means]
    # The paired lower is independent of the nineteen-site angles.
    largeopts=[options(v,False) for v in ys[867][:2]]
    if any(not l for l in smallopts+largeopts):
        return None
    a=min((loss(17,p,q)+loss(26,z,q)+loss(29,z,p)).lower() for p,q,z in product(*smallopts))
    b=min(loss(867,i,j).lower() for i,j in product(*largeopts))
    return a+b+hcost


def main():
    # Any other zero costs more than this even before reference quantization.
    threshold=arb(".00215")
    cap=int((threshold/lam).sqrt().upper().ceil().unique_fmpz())
    assert lam*cap**2>threshold
    count=arb(cap).zeta_nzeros().unique_fmpz()
    assert count is not None
    for n,z in enumerate(acb.zeta_zeros(1,int(count)),1):
        if n not in relaxed:
            assert m["paired_lower"](m["state"](z.imag))+lam*z.imag**2>threshold
    # At a continuous minimum the first variation vanishes. Bounding the second
    # derivative shell-by-shell gives a total excess <=88*(pi/(2m))^2.
    tail=2300
    assert relaxed[255]+88*(pi/(2*tail))**2<relaxed[25]<threshold
    hard=[]
    for k in range(343,tail):
        u=trial(255,k,False)
        if u>=relaxed[25]:
            u=trial(255,k,True)
        assert u<threshold
        if u<relaxed[25]:
            continue
        n=255 if k>343 else 25
        hi=u if n==255 else trial(25,k,True)
        competitor=25 if n==255 else 255
        lo=lower(competitor,k,hi)
        assert lo is None or lo>hi
        hard.append((k,n))
    assert hard==[(343,25),(373,255),(374,255)]
    print("finite orders needing a quantized competitor comparison",hard)
    print("last exceptional order",343,"rigorous uniform quadratic tail from",tail)


if __name__=="__main__":
    main()
