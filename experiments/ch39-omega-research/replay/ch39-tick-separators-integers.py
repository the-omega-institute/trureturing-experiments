"""Replay the lower bounds using only integers and fractions, given certified input cells."""
from fractions import Fraction as Q
from math import comb, factorial, lcm
from pathlib import Path
from runpy import run_path
import json


def chebyshev_upper(base, ticks):
    den=1<<72
    low=high=den
    x,y=base,base+1
    a,b=x,y
    for _ in range(1,ticks):
        corners=[r*s for r in (x,y) for s in (a,b)]
        lower=2*min(corners)//den-high
        upper=-((-2*max(corners))//den)-low
        low,high=a,b
        a,b=lower,upper
    return a,b


def mul_linear(poly,constant,slope):
    out=[0]*(len(poly)+1)
    for i,a in enumerate(poly):
        out[i]+=a*constant
        out[i+1]+=a*slope
    return out


def bernstein_positive(coefficients, depth=16):
    if min(coefficients)>0:
        return
    assert depth>0
    left,right=[],[]
    row=coefficients
    degree=len(row)-1
    while row:
        scale=1<<(len(row)-1)
        left.append(row[0]*scale)
        right.append(row[-1]*scale)
        row=[a+b for a,b in zip(row,row[1:])]
    bernstein_positive(left,depth-1)
    bernstein_positive(right[::-1],depth-1)


def tail_positive(a, peak, ceiling, h):
    deg=len(a)-1
    u,v=1<<800,h*h
    D=[1]
    slopes=[j*j*u for j in range(1,deg+1)]
    for m in slopes:
        D=mul_linear(D,v,m)
    leading=(ceiling-peak)*(1<<40)-a[0]
    scale=leading.denominator
    poly=[leading.numerator*x for x in D]
    for index,m in enumerate(slopes,1):
        # D/(v+m*x), by exact division from the constant term up.
        last=0
        for k in range(deg):
            value,rest=divmod(D[k]-last*m,v)
            assert rest==0
            poly[k]-=scale*a[index]*v*value
            last=value
        assert last*m==D[-1]
    # Multiply every Bernstein coefficient by deg! to keep integers.
    factors=[factorial(j)*factorial(deg-j) for j in range(deg+1)]
    bern=[sum(poly[j]*comb(i,j)*factors[j] for j in range(i+1)) for i in range(deg+1)]
    bernstein_positive(bern)


def torus_upper(c, tick):
    special={1:Q(".369381"),5:Q(".436809"),7:Q(".925928")}
    if tick not in special:
        return sum(map(abs,c))
    peak=special[tick]
    a,b,d,e=c
    assert peak>abs(a)
    # Coefficients of (peak-a*x)^2 - (b+(d+e)*x)^2 -(d-e)^2*(1-x^2).
    p0=peak**2-b*b-(d-e)**2
    p1=-2*(peak*a+b*(d+e))
    p2=a*a-(d+e)**2+(d-e)**2
    assert min(p0-p1+p2,p0+p1+p2)>0
    if p2>0 and -2*p2 < -p1 < 2*p2:
        assert 4*p2*p0>p1*p1
    return peak


def certify(tick, rec, prefix, levels, bits):
    a=rec["a"]; f=1<<40
    c=[Q(v,1<<rec["depth"]) for v in rec["c"]]
    ceiling=Q(".0001" if tick==5 else ".002")
    threshold=ceiling+sum(map(abs,a))/Q(f*2**bits)+sum(map(abs,c))*Q(1,2**16)
    # Each target is in [l,l+1]/2^288.
    target=sum(t*(l+int(t<0)) for t,l in zip(a,levels))/Q(f*2**288)
    assert target>threshold
    cap=int(ceiling*(1<<112))
    native_pad=sum(map(abs,c))*Q(1,2**96)
    extra=-((-native_pad*(1<<112))//1)
    cs=rec["c"]
    for h,*cosbins in prefix:
        high=0
        for j,t in enumerate(a):
            h2=(h+int(t>=0))**2
            high-=((-t*h2*(1<<72))//(h2+(j*j<<800)))
        for t,x in zip(cs,cosbins):
            lo,hi=chebyshev_upper(x,tick)
            high-=((-t*(hi if t>=0 else lo)*(1<<40))//(1<<rec["depth"]))
        assert high+extra<cap
    tail_positive(a,torus_upper(c,tick),ceiling,prefix[-1][0])
    print("tick",tick,"excluded at time precision",bits,"offsets 0..",len(a)-1)


def main():
    root=Path(__file__).resolve().parent
    cells=json.loads((root/"ch39-tick-separators.json").read_text())
    source=json.loads((root/"ch39-minimum-tick-finite-inputs.json").read_text())
    horizon=run_path(str(root/"ch39-prime-clock-height-horizon.py"))["main"]
    p=cells["prefix"]
    assert len(p)==3200 and all(x[0]+1<y[0] for x,y in zip(p,p[1:]))
    assert p[-1][0]+1<2**415
    for tick in range(1,8):
        horizon(tick,15)
        certify(tick,cells["dual"][str(tick)],p,source["levels"],96 if tick==1 else 272)
    run_path(str(root/"ch39-minimum-tick-finite-replay.py"))["main"]()
    print("Minimum common positive shell-clock tick on the first 54 readings: two at 96 bits and eight at 272 bits.")


if __name__=="__main__":
    main()

