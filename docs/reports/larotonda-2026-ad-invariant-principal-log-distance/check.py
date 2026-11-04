#!/usr/bin/env python3
"""Exact checks for the conjecture following Larotonda 2610.31535v1, Prop. 4.3.
Angles are rational multiples of pi; squared lengths are reported in units pi^2.
This checks algebra and explicit curves, not a Lean proof or an exact distance.
"""
from fractions import Fraction as Q
from itertools import product
import sympy as s


def norm_sq(a, c):
    return sum(x*x for x in a) + c*sum(a)**2


def shift(a, j=0):
    b = list(a)
    b[j] -= 2
    return tuple(b)


def same_exp(a, b):
    assert all(((x-y)/2).denominator == 1 for x, y in zip(a, b))
    assert all(s.simplify(s.exp(s.I*s.pi*s.Rational(x.numerator,x.denominator)) -
                          s.exp(s.I*s.pi*s.Rational(y.numerator,y.denominator))) == 0
               for x, y in zip(a, b))


# Fidelity: Section 4's scalar SU(5) example exp(2*pi*i/5) I.
a5 = (Q(2,5),)*5
assert sum(a5) == 2 and max(map(abs, a5)) == Q(2,5)
b5 = shift(a5)
same_exp(a5, b5)
assert sum(b5) == 0 and max(map(abs, b5)) == Q(8,5)
# A traceless logarithm has eigenangles (2/5 + 2*k_j)*pi and sum k_j=-1.
# Consequently some k_j <= -1, and its absolute eigenangle is at least 8/5*pi.
for k in range(-25, 0):
    assert abs(Q(2,5)+2*k) >= Q(8,5)
print('FIDELITY Section 4 SU(5): principal op/pi=2/5, trace/(i*pi)=2;')
print('  traceless shifted log op/pi=8/5, trace=0; same exponential checked exactly')

# Fidelity: Theorem 2.1 for diagonal Frobenius norm (c=0), exact small cases.
fidelity_cases = 0
for a in product([Q(-1), Q(-1,2), Q(0), Q(1,2), Q(1)], repeat=2):
    for k in product(range(-3,4), repeat=2):
        b = tuple(x+2*j for x,j in zip(a,k))
        assert norm_sq(a,Q(0)) <= norm_sq(b,Q(0))
        fidelity_cases += 1
print(f'FIDELITY Theorem 2.1: {fidelity_cases} diagonal Frobenius comparisons passed')

# The required n=2 strict counterexample, c=4.
a2 = (Q(4,5), Q(1,2)); b2 = shift(a2); c = Q(4)
same_exp(a2, b2)
assert max(map(abs,a2)) == Q(4,5) < 1
assert norm_sq(a2,c) == Q(153,20)
assert norm_sq(b2,c) == Q(73,20)
assert norm_sq(b2,c) < norm_sq(a2,c)
print('REFUTATION n=2 c=4: H/pi=(4/5,1/2), Hprime/pi=(-6/5,1/2)')
print('  principal length^2/pi^2=153/20; competitor length^2/pi^2=73/20; gap=-4')

# Cross-check the search seat's n=3 data and its whole c>3 threshold.
c_sym = s.symbols('c', real=True)
a3 = tuple(map(Q, ['2/5','9/20','7/20'])); b3 = shift(a3)
same_exp(a3,b3)
assert s.simplify(norm_sq(b3,c_sym)-norm_sq(a3,c_sym) + 4*(c_sym-3)/5) == 0
assert norm_sq(a3,Q(4)) == Q(1249,200)
assert norm_sq(b3,Q(4)) == Q(1089,200)
assert max(map(abs,a3)) == Q(9,20) < Q(1,2)
print('CROSSCHECK n=3: threshold c>3; c=4 gives principal 1249/200, competitor 1089/200')

# General rank-one winding-shift identity (theta/pi coordinates).
x,y = s.symbols('x y', real=True)
assert s.expand(norm_sq((x-2,y),c_sym)-norm_sq((x,y),c_sym) -
                4*((1-x)+c_sym*(1-x-y))) == 0
# Positive and negative c branches give witnesses for every admissible nonzero c.
eps_pos = c_sym/(4*(1+c_sym))
eps_neg = -c_sym/(4*(1+c_sym))
for eps, yval, gap in [(eps_pos,s.Rational(1,2),-c_sym),
                       (eps_neg,-s.Rational(1,2),c_sym)]:
    aa=(1-eps,yval); bb=(aa[0]-2,aa[1])
    assert s.factor(norm_sq(bb,c_sym)-norm_sq(aa,c_sym)-gap) == 0
print('FAMILY symbolic: c>0 gap/pi^2=-c; -1/n<c<0 gap/pi^2=c, n>=2')

# Positive definiteness and normalization: Q_c(H) splits into traceless + central.
h11,h22,hr,hi = s.symbols('h11 h22 hr hi',real=True)
H=s.Matrix([[h11,hr+s.I*hi],[hr-s.I*hi,h22]])
T=s.trace(H); H0=H-(T/2)*s.eye(2)
qH=s.trace(H*H)+c_sym*T*T
assert s.expand(qH-(s.trace(H0*H0)+(s.Rational(1,2)+c_sym)*T*T)) == 0
# The coefficient matrix for the real diagonal part is I+c J, eigenvalues 1,1+2c.
G=s.eye(2)+c_sym*s.ones(2)
assert set(G.eigenvals()) == {s.Integer(1),1+2*c_sym}
# Nontrivial exact unitary conjugation, testing trace and square-trace invariance.
U=s.Matrix([[1,s.I],[s.I,1]])/s.sqrt(2)
assert s.simplify(U.conjugate().T*U) == s.eye(2)
HU=s.simplify(U*H*U.conjugate().T)
assert s.expand(s.trace(HU)-T) == 0
assert s.simplify(s.trace(HU*HU)-s.trace(H*H)) == 0
print('ADMISSIBILITY n=2: quadratic split and nontrivial unitary invariance checked symbolically')
print('  c=4: Q=Tr(H^2)+4(Tr H)^2 >= Tr(H^2) >= op(H)^2')

family_count=0
for n in range(2,9):
    for cc in [Q(-1,n)+Q(1,100*n), Q(-1,2*n), Q(-1,100*n),
               Q(1,100),Q(1,2),Q(1),Q(4),Q(100)]:
        assert cc > Q(-1,n) and cc != 0 and 1+cc > 0
        eps=abs(cc)/(4*(1+cc))
        aa=(1-eps,Q(1,2) if cc>0 else Q(-1,2))+(Q(0),)*(n-2)
        bb=shift(aa)
        assert 0<eps<Q(1,4)
        assert max(map(abs,aa)) < 1
        assert norm_sq(bb,cc)-norm_sq(aa,cc) == -abs(cc) < 0
        # For c<0 multiply F by 1/sqrt(1+n*c); for c>0 no rescaling.
        scale_sq=1/min(Q(1),1+n*cc)
        assert scale_sq > 0
        assert scale_sq*norm_sq(bb,cc) < scale_sq*norm_sq(aa,cc)
        family_count+=1
print(f'FAMILY: {family_count} exact parameter tests, dimensions 2..8, both signs of c')

# Explicit shorter curve: gamma(t)=diag(exp(i*pi*t*b_j)).
# Check endpoints and the literal left logarithmic derivative gamma^* gamma'.
t=s.symbols('t',real=True)
gamma=s.diag(*[s.exp(s.I*s.pi*t*s.Rational(v.numerator,v.denominator)) for v in b2])
assert s.simplify(gamma.subs(t,0)) == s.eye(2)
endpoint=s.diag(*[s.exp(s.I*s.pi*s.Rational(v.numerator,v.denominator)) for v in a2])
assert all(s.simplify(v)==0 for v in gamma.subs(t,1)-endpoint)
velocity=s.simplify(gamma.conjugate().T*gamma.diff(t))
assert velocity == s.diag(*[s.I*s.pi*s.Rational(v.numerator,v.denominator) for v in b2])
speed_sq=s.Rational(73,20)*s.pi**2
assert s.integrate(s.sqrt(speed_sq),(t,0,1)) == s.sqrt(speed_sq)
print('CURVE: unitary, endpoints I and exp(iH); left velocity=iHprime;')
print('  length=pi*sqrt(73/20) < pi*sqrt(153/20), entirely exact')
print('ALL_CHECKS_OK')
