#!/usr/bin/env python3
"""Exact checks of Hagihara arXiv:2611.10920v1, Conjecture 3.3.

Run with Python >= 3.12 and SymPy. No floating-point arithmetic is used.
Finite recurrence checks are evidence for the displayed formulas, not an
infinite positivity proof or a Lean kernel certificate.
"""
from fractions import Fraction as F
from itertools import product
from math import comb
import hashlib
import sys
import sympy as s

assert sys.version_info >= (3, 12)
x, z, b, y = s.symbols('x z b y', real=True)
r = s.symbols('r', real=True)
D = (3*x+11)*(3*x+8)
P = -(3*x*x+19*x+29+r)/D
Q = (x*x+4*x+4+r)*(x+4)/(D*(3*x+5))
T = -s.Rational(1,3)-1/(9*x*x)

# Fidelity: source Example 3.2, the first-order shifted sum f1 and its
# displayed second-order recurrence (source TeX label ex:2f1<f2).
def franel(n):
    return 2*n*sum(comb(n,k)**3 for k in range(n+1))
def U(n):
    return F(7*n*n+21*n+16,(n+1)*(n+2))
def V(n):
    return F(8*(n+1)**2,n*(n+2))
for n in range(1,31):
    assert franel(n+2) == U(n)*franel(n+1)+V(n)*franel(n)
print('FIDELITY: Example 3.2 Franel recurrence, n=1..30: OK')

residual = ((x+z)*(x+z+1)*(x+z+2)
    -(3*x+11)*(x+z)*(x+z+1)
    +(3*x*x+19*x+29+b*b)*(x+z)
    -(x*x+4*x+4+b*b)*(x+4))
assert s.expand(residual-(z-4)*((z-2)**2+b*b)) == 0
assert s.simplify(s.limit(P,x,s.oo)+s.Rational(1,3)) == 0
assert s.simplify(s.limit(x*x*(P+s.Rational(1,3)),x,s.oo)
                  -(s.Rational(1,27)-r/9)) == 0
assert s.expand((x+3)**2-((x+2)**2+r)) == 2*x+5-r
print('FAMILY: polynomial residual and factorwise norm certificate 2*n+5-r: OK')

M = s.Matrix([[1,1,0],[s.Rational(4,5),s.Rational(2,5),b/5],
              [s.Rational(1,2),(6-b*b)/40,b/8]])
assert s.factor(M.det()) == b*(b*b+4)/200
M2=M.subs(b,2)
expected_inverse=s.Matrix([[1,-s.Rational(25,8),5],
                           [0,s.Rational(25,8),-5],
                           [-2,s.Rational(45,8),-5]])
assert M2.inv() == expected_inverse
assert M2.det() == s.Rational(2,25)
assert s.Rational(1)-s.Rational(239,512) == s.Rational(273,512)
print('BASIS: det(M_b)=b*(b^2+4)/200; b=2 inverse and 273/512 margin: OK')

gap=((9*r-12)*x*x-57*x-88)/(9*x*x*D)
assert s.cancel(T-P-gap) == 0
assert s.limit(x*x*gap,x,s.oo) == (3*r-4)/27
cert=s.cancel((gap.subs(r,4)-1/(16*x*x))*(144*x*x*D))
assert s.expand(cert) == 303*x*x-1425*x-2200
assert s.expand(cert.subs(x,y+10)) == 303*y*y+4635*y+13850
assert s.expand(z**3-z**2+z/3-s.Rational(1,27)-(z-s.Rational(1,3))**3)==0
print('BOUNDARY: exact T-P formula, positive limit for r>4/3, x>=10 certificate: OK')

# Exact small cases for all three real basis solutions (complex rational
# arithmetic at b=2); the recurrence is evaluated independently.
def pq(n,rr):
    return (F(-(3*n*n+19*n+29+rr),(3*n+11)*(3*n+8)),
            F((n*n+4*n+4+rr)*(n+4),(3*n+11)*(3*n+8)*(3*n+5)))
t=[F(1)]; gr=[F(1)]; gi=[F(0)]
for n in range(17):
    t.append(t[-1]*F(n+4,3*n+5))
    gr.append((gr[-1]*(n+2)-gi[-1]*2)/F(n+4))
    gi.append((gi[-1]*(n+2)+gr[-2]*2)/F(n+4))
for f in [t,[a*c for a,c in zip(t,gr)],[a*c for a,c in zip(t,gi)]]:
    for n in range(15):
        p,q=pq(n,4)
        assert f[n+3] == f[n+2]+p*f[n+1]+q*f[n]
assert t[:5] == [F(1),F(4,5),F(1,2),F(3,11),F(3,22)]
assert s.simplify(s.prod(1-(2-2*s.I)/k for k in range(1,4))+5*s.I/3)==0
print('PRODUCTS: three b=2 basis solutions, n=0..14; baseline and Pochhammer denominator: OK')

# Direct forward recurrence in exact rationals for b^2=2. This uses no
# closed-form product solutions. Two arbitrary real perturbation patterns
# are represented by rational initial values within 1/1024 of the baseline.
# Fractions are reduced exactly after every operation; three values suffice.
for name, signs in [('plus-minus-plus',(1,-1,1)),('minus-plus-minus',(-1,1,-1))]:
    init=[F(1),F(4,5),F(1,2)]
    f0,f1,f2=[a+F(v,2048) for a,v in zip(init,signs)]
    assert min(f0,f1,f2)>0
    digest=hashlib.sha256()
    for n in range(9998):
        p,q=pq(n,2)
        fn=f2+p*f1+q*f0
        assert fn>0, (name,n+3)
        f0,f1,f2=f1,f2,fn
        # Byte serialization avoids decimal conversion limits for large ints.
        for a in [fn.numerator,fn.denominator]:
            raw=a.to_bytes((a.bit_length()+7)//8,'big')
            digest.update(len(raw).to_bytes(8,'big'));digest.update(raw)
    print(f'DIRECT_RECURRENCE: r=2, {name}, delta=1/2048, n=0..10000 positive; '
          f'last numerator bits={f2.numerator.bit_length()}, denominator bits={f2.denominator.bit_length()}; '
          f'trace_sha256={digest.hexdigest()}',flush=True)
print('ALL_CHECKS_OK')
