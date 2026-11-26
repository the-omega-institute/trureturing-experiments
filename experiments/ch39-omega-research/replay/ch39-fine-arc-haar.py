"""Independent angles cut out of one original arc, with exact rational shell coordinates."""
from fractions import Fraction
from math import isqrt, lcm
from flint import arb, ctx, fmpz

ctx.prec=160
theta=(1-arb(12)/1734).acos()
h=Fraction(1,100)
n=12
K=int((4*arb.pi()/theta*max(n,arb(h.denominator)/h.numerator)).ceil().unique_fmpz())
cap=theta*arb(h.numerator)/h.denominator
primes=[]
p=17
while len(primes)<n:
    p+=2
    if p%8==1 and fmpz(p).is_prime():primes.append(p)


def multiply(z,w):
    a,b=z
    c,d=w
    return a*c-2*b*d,a*d+b*c


def power(z,k):
    out=(1,0)
    while k:
        if k%2:out=multiply(out,z)
        z=multiply(z,z)
        k//=2
    return out


point=(Fraction(25),Fraction(11))
D=1
total=arb(0)
upper_bits=0
assert 4*arb.pi()/K <= cap
for p in primes:
    a,b=next((isqrt(p-2*b*b),b) for b in range(1,isqrt(p//2)+1)
             if isqrt(p-2*b*b)**2==p-2*b*b)
    rotation=multiply((a,b),(a,b))
    angle=2*(arb(2).sqrt()*b/a).atan()
    for k in range(1,K+1):
        t=k*angle/(2*arb.pi())
        frac=t-t.floor()
        arg=2*arb.pi()*(frac if frac<arb(".5") else frac-1)
        if abs(arg)<2*arb.pi()/K and not arg.contains(0):break
    else:raise AssertionError("Dirichlet turn bound")
    num=power(rotation,k)
    rot=(Fraction(num[0],p**k),Fraction(num[1] if arg>0 else -num[1],p**k))
    assert rot[0]**2+2*rot[1]**2==1
    point=multiply(point,rot)
    assert point[0]**2+2*point[1]**2==867
    D=lcm(D,point[0].denominator,point[1].denominator)
    total+=abs(arg)
    upper_bits+=k*p.bit_length()
assert 0<total and total<theta/2 and D.bit_length()<=upper_bits+1
print("prime clocks",primes)
print("Dirichlet bound K",K,"total original arc fraction",(total/theta).upper().str(12))
print("integer dilation bit length for the selected vertices",D.bit_length())
