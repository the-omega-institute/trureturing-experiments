"""Sharp singular-value hierarchy for a real linear clock plus the other clock squared."""
from fractions import Fraction as F
from flint import arb
from sympy import symbols,factor,cancel

t,A=symbols("t A")
s2=8*A**2/(1+4*A**2)
d=s2*t*(1-t)**2/(1-s2*t)
for sign in (1,-1):
    test=d+sign*A*(2*t-1)-A**2
    polynomial=4*A**3+sign*(4*A**2-8*A**2*t)+A*(4*t**2-6*t+1)+sign
    assert cancel(test-A*(2*t-1)*polynomial/(1+4*A**2-8*A**2*t))==0
for c2 in (F(81,98),F(8741547,12563572)):
    c2=arb(c2.numerator)/c2.denominator
    threshold=((1-c2)/(1+c2)).sqrt()/2
    assert threshold>arb(".01947") and threshold<arb(1)/4
    print("minimum second/first singular value",threshold.str(15),
          "attaining linear-to-quadratic coefficient ratio",(2/(1+c2)).sqrt().str(15))
