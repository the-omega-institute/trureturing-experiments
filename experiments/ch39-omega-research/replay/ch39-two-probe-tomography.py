"""Complete exceptional set for tomography of the fixed 50-site register."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
s=run_path(str(root/"ch39-fixed-scale-common-reading.py"))
B,base,rounded,multiply=s["B"],s["base"],s["rounded"],s["multiply"]


def times_linear(poly,x):
    result=[B(0) for _ in range(len(poly)+1)]
    for i,v in enumerate(poly):
        result[i]=result[i]+rounded(v*x)
        result[i+1]=result[i+1]+v
    return list(map(rounded,result))


def extra_polynomial(tri,c,models):
    # c*(T_plus-T_minus)(eta) = W(s)*P(s)/prod(s+r_i), s=eta/c.
    # Interpolate P(u*s0) using five negative poles, with u scaled for clarity.
    unit=10**6
    offsets=[Q(j*j) for j in base["OFFSETS"] if j]+[s["eta"]/c]
    terms=[]
    r=tri.poles
    diff={i:(b[1]-a[1]) for i,a,b in zip(tri.ids,models[0].ref,models[1].ref)}
    nodes=(45,81,127,210,473)
    for i in nodes:
        v=r[i]
        numerator=multiply(diff[i],*(r[k]-v for k in tri.ids if k!=i))
        # The rational function uses W(s)=prod(s-x); at -v all 44 signs cancel.
        den=multiply(*(v+B(x) for x in offsets),*((r[k]-v)/B(unit) for k in nodes if k!=i))
        poly=[rounded(numerator/den)]
        for k in nodes:
            if k!=i:poly=times_linear(poly,r[k]/B(unit))
        terms.append(poly)
    return [rounded(sum(row,B(0))) for row in zip(*terms)]


def main():
    tri=s["ThreeFits"]()
    c=Q(base["data"]("ch39-fixed-scale-key.json")[0]["c"])
    _,models=tri.models(c)
    poly=extra_polynomial(tri,c,models)
    a,b,d,e,f=poly
    assert a.lo>0 and b.hi<0 and min(d.lo,e.lo,f.lo)>0
    assert b.lo*b.lo<4*a.lo*d.lo
    # The first three coefficients form a positive definite quadratic. The
    # two highest terms are positive on the entire positive half-line.
    offsets={c*j*j for j in base["OFFSETS"] if j}
    assert len(offsets)==43 and s["eta"] not in offsets
    print("Known 50 heights and scale: two probes uniquely recover the weights except at 44 positive second spacings.")
    scale=max(max(abs(v.lo),abs(v.hi)) for v in poly)
    print("Normalized quartic coefficients:",[(float(v.lo/scale),float(v.hi/scale)) for v in poly])


if __name__=="__main__":
    main()
