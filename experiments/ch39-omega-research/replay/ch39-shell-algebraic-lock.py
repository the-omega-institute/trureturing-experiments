"""Height-denominator obstruction when three algebraic phase locks are imposed on the two shells."""
from fractions import Fraction as F
from flint import arb,acb,ctx

# At one prime above 17 in Q(sqrt(-42),sqrt(-93)) the two roots reduce to 3 and -3.
assert (3**2+42)%17==0 and ((-3)**2+93)%17==0
assert (11+2*3)%17==0 and (11-2*3)%17!=0 and (14-3)%17!=0
assert 11**2+4*42==17**2==14**2+93
# General supplied field degree up to 10^50 for linear three-deficit selectors.
ctx.prec=540
h=acb.zeta_zero(719).imag
B=10**140
z=(h*B).floor().unique_fmpz()
assert z is not None
center=F(int(z),B)
N=16*10**50
def farey_pair(x,maximum):
    n,d=x.as_integer_ratio()
    p0,q0,p1,q1=0,1,1,0
    while d:
        a=n//d
        if q0+a*q1>maximum:
            break
        p0,q0,p1,q1=p1,q1,p0+a*p1,q0+a*q1
        n,d=d,n-a*d
    k=(maximum-q0)//q1
    return sorted((F(p0+k*p1,q0+k*q1),F(p1,q1)))
lo,hi=farey_pair(center,N)
assert hi.numerator*lo.denominator-lo.numerator*hi.denominator==1
assert max(lo.denominator,hi.denominator)<=N<lo.denominator+hi.denominator
assert arb(lo.numerator)/lo.denominator<h<arb(hi.numerator)/hi.denominator
print("gamma719 Farey neighbor denominators and minimum interior denominator",
      lo.denominator,hi.denominator,lo.denominator+hi.denominator)
print("prime valuation of exp(2i omega_17,1)",4,
      "excluded coefficient field degree up to",10**50)
