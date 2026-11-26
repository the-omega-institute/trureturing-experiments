"""Quadratic characters for eight independent clocks and constants of small-loss asymptotics."""
from fractions import Fraction
from math import isqrt
from pathlib import Path
from runpy import run_path
from flint import arb,acb

model=run_path(str(Path(__file__).with_name("ch39-shell-third-shell.py")))
chars=[]
volumes=[]
for r,(ds,c,_) in model["shells"].items():
    for d in ds:
        cosine=Fraction(r-d,r)
        assert (2*cosine).denominator>1
        a,b=cosine.as_integer_ratio()
        rad=b*b-a*a
        k=max(i for i in range(1,isqrt(rad)+1) if rad%(i*i)==0)
        chars.append(-rad//(k*k))
    def integrand(x,analytic):
        s=(2*x).sin()
        h=2-2*c*s+(c*c-arb(1)/2)*s*s
        return 1/h.sqrt(analytic=analytic)
    v=acb.integral(integrand,0,2*arb.pi())/2
    assert v.imag.contains(0) and v.real>0 and v.real.rad()<arb("1e-16")
    volumes.append(v.real)
assert len(set(chars))==12
print("distinct quadratic characters",chars)
pi=arb.pi()
constant=arb(1)/(32*pi**6)
for v in volumes:
    constant*=v
print("four quartic sublevel areas",*[v.str(14) for v in volumes])
print("fixed odd orders, coefficient of epsilon^2",constant.str(16))
print("continuous phases, elementary lower coefficient", (4/(arb(3)**2*arb(3).sqrt()*pi**2)).str(15))
print("paired-loss upper asymptotic coefficient",(8*arb(2).sqrt()/(3*pi**5)).str(15))
