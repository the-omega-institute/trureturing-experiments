"""Enclose the full mixing-modulus fiber of the measured two-matrix spectral pencil."""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_poly

folder=Path(__file__).resolve().parent
m=run_path(str(folder/"ch39-shell-gw-mixed-determinant-cp.py"))
density=run_path(str(folder/"ch39-shell-gw-positive-density-cp.py"))


def invisible(values):
    a,b,c=[v**2 for v in values]
    return [(c-b)/(c-a),arb(-1),(b-a)/(c-a)]


def fiber():
    V=m["ckm"](arb("1.139"))
    r,s=invisible(m["up"]),invisible(m["down"])
    return [[arb_poly([abs(V[i,j])**2,r[i]*s[j]]) for j in range(3)] for i in range(3)]


def isolate(poly,lo,hi):
    lo,hi=arb(lo),arb(hi)
    sign=1 if poly(hi)>0 else -1
    assert sign*poly(lo)<0 and sign*poly(hi)>0
    for _ in range(64):
        mid=(lo+hi)/2
        if sign*poly(mid)>0:
            hi=mid
        else:
            assert sign*poly(mid)<0
            lo=mid
    return lo.union(hi)


def angles(B,t):
    s13=B[0][2](t).sqrt()
    denom=1-s13**2
    return [(B[0][1](t)/denom).sqrt().asin(),s13.asin(),(B[1][2](t)/denom).sqrt().asin()]


def target_bound():
    B=fiber()
    a,b,c=[B[0][i]*B[1][i] for i in range(3)]
    polynomial=4*a*b-(c-a-b)**2
    # The quartic coefficient vanishes identically for this rank-one direction:
    # 4 s0² s1² = (s2²-s0²-s1²)² since s0+s1+s2=0.
    assert polynomial[4].contains(0) and polynomial[3]>0
    lo=isolate(polynomial,"-.020","-.019")
    hi=isolate(polynomial,".044",".045")
    assert polynomial(arb(251))<0 and polynomial(arb(253))>0
    for row in B:
        for p in row:
            assert p(lo)>0 and p(hi)>0
    assert B[1][0](arb(251))<0
    assert polynomial(arb(0))>0
    print("physical endpoints",lo.str(11),hi.str(11))
    for t in (lo,hi):
        th=angles(B,t)
        t12,t13,t23=th
        odd=2*t12.sin()*t12.cos()*t23.sin()*t23.cos()*t13.sin()
        base=(t12.sin()*t23.cos())**2+(t12.cos()*t23.sin()*t13.sin())**2
        print("angles",*[v.str(8) for v in (t12,t23,t13)],"cos phase",((B[1][0](t)-base)/odd).str(7))
    # Monotone ratios of affine functions attain their extrema at the endpoints.
    limits=[a.union(b) for a,b in zip(angles(B,lo),angles(B,hi))]
    t12,t13,t23=m["angles"]
    d,s,b=m["down"]
    hs=(d*d+s*s+b*b).sqrt()
    angular12=abs(limits[0]-t12).upper()
    angular23=abs(limits[2]-t23).upper()
    sine=limits[1].sin().upper()+t13.sin()
    cos_error=abs(limits[1].cos()-t13.cos()).upper()
    error=angular12*(d*d+s*s).sqrt()+hs*(angular23+(sine*sine+cos_error*cos_error).sqrt())
    print("joint target change below",error.upper().str(9))
    return error


if __name__=="__main__":
    error=target_bound()
    for n,width in ((43,"3e-8"),(195,"1e-8")):
        op=density["window"](n,arb(width))
        density["cp"]["check"](n,op,error)
