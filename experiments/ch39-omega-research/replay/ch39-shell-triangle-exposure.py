"""Density of exposed Wilson readout triples on the two arithmetic shells."""
from fractions import Fraction as F
from flint import acb,arb,ctx
ctx.prec=112
miss=[]
for c2 in (F(81,98),F(8741547,12563572)):
    c=(arb(c2.numerator)/c2.denominator).sqrt()
    u=-c.log()
    fail=(4/arb.pi()**2*acb.integral(lambda z,a:1/(acb(0,1)*z).sinc(),0,u)).real
    assert 0<fail<1
    miss.append(fail)
    print("axis cosine and exposed fraction",c.str(14),(1-fail).str(15))
print("at least one and both shells individually expose the target",
      (1-miss[0]*miss[1]).str(15),((1-miss[0])*(1-miss[1])).str(15))
