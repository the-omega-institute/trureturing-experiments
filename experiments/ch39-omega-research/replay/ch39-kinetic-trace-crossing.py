"""Trace balance and required complementary scalar exposure on the original 18-mode window."""
from fractions import Fraction
from pathlib import Path
import runpy

from flint import arb

root=Path(__file__).resolve().parent
model=runpy.run_path(str(root/"ch39-quaternion-gauge.py"))
q=model["q"]
r=q("31.0740").union(q("31.6548"))
d=q(460)**2/q("1.04")
ys=[x/r**2 for x in model["xs"]]
low_trace=sum((-x).exp() for x in ys)
low_first=sum(x*(-x).exp() for x in ys)
# Continuum original inertia, geometric zero-35 square. Above the 18-low
# trial space all energies exceed 100^2 by the Schur bound below first pole.
k=100*arb.pi()/d.sqrt()
assert 959990-10000-36*k/k.sin()>0
first=10000/r.upper()**2
c=d/(3*r.upper()**2)
# Both uniformly summable tails use the boundary-zero subspace (codim 38).
# For the first moment, y exp(-y) <= exp(-y/2) once y>=2.
assert c>2 and first>1
x=(-c).exp()
z=(-c/2).exp()
trace_upper=low_trace.upper()+20*(-first).exp()+270*x/(1-x)
moment_upper=low_first.upper()+20*first*(-first).exp()+270*z/(1-z)
mu=moment_upper/low_trace.lower()
light_weight=(1-mu)*low_trace.lower()
exposure=-low_trace.lower()*mu.log()/2
assert 0<mu and mu<1
# A positive heat trace of R states matched in value and slope has the same
# mean y. Counting raw logarithmic exposures bounds them from below.
for name,v in (("left trace lower",low_trace.lower()),("left trace upper",trace_upper),
               ("mean y upper",mu),("required sub-Q trace lower",light_weight),
               ("right log exposure lower",exposure),
               ("weak minus color one-loop shift lower",exposure/(6*arb.pi())),
               ("color minus hypercharge shift lower",exposure/(10*arb.pi()))):
    print(name,v.lower().str(12),v.upper().str(12))
# Exact scalar trace indices for (left doublets, right singlets, Higgs, triplet).
L,R,H,C=[tuple(map(Fraction,row)) for row in (
    ("2/5",2,1),("8/5",0,1),("3/10","1/2",0),("1/5",0,"1/2"))]
assert [a+b for a,b in zip(L,R)]==[2]*3 and [a+b for a,b in zip(H,C)]==[Fraction(1,2)]*3
assert ys[0].upper()+q(2).log()<1 and q(3).log()>1
