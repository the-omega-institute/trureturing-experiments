#!/usr/bin/env python3
"""Exact checks for the Debus–Schabert uniform Reznick counterexample.

Python 3 plus SymPy; every mathematical comparison uses integers or rationals.
This program does not check the infinite-dimensional analytic argument or a Lean proof.
"""
from fractions import Fraction as Q
from itertools import combinations, product
import json
import sympy as s

DEN = 10**9
PARTITIONS = [(9,0,0),(8,1,0),(7,2,0),(7,1,1),(6,3,0),(6,2,1),
              (5,4,0),(5,3,1),(5,2,2),(4,4,1),(4,3,2),(3,3,3)]
NUMERATORS = [1000000000,234635012,178775790,7016415,166340945,5047384,
              163611602,4928454,3677868,4910485,3646976,3615007]
MU = dict(zip(PARTITIONS, (Q(n,DEN) for n in NUMERATORS)))
EXPECTED_FLOORS = [
    [234635012,4837568,83691,86419,187694598,29913962,879731,39334,
     26850362,470290,1957,330440,158345,158016,69021606],
    [7016415,1297527,30401,1526629,1189725,159087,202848,11454,4701,20630]]


def p(xs, j):
    return sum((x**j for x in xs), s.Integer(0))


def rob(a,b,c):
    return (a**6+b**6+c**6-a**4*b**2-a**4*c**2-b**4*a**2
            -b**4*c**2-c**4*a**2-c**4*b**2+3*a**2*b**2*c**2)


def hp(xs):
    return (p(xs,10)*p(xs,4)**2+p(xs,10)*p(xs,8)+p(xs,6)**3
            -p(xs,12)*p(xs,6)-2*p(xs,8)*p(xs,6)*p(xs,4))/2


def ht(xs):
    return sum(((a*b*c)**4*rob(a,b,c) for a,b,c in combinations(xs,3)),s.Integer(0))


def ident(a,b):
    assert s.expand(a-b) == 0


# Fidelity: the power-sum identity in source Example 6.6, before summing over i.
P2,P4,P6,P8,z=s.symbols('P2 P4 P6 P8 z')
def A(a,b,c,d):
    return a**5/s.Integer(18)+3*d*a+6*c*b-3*c*a**2
shifted = A(P2-z**2,P4-z**4,P6-z**6,P8-z**8)
stated = (A(P2,P4,P6,P8)+(-5*P2**4/s.Integer(18)+6*P2*P6-3*P8)*z**2
          +(5*P2**3/s.Integer(9)-9*P6)*z**4
          +(22*P2**2/s.Integer(9)-6*P4)*z**6
          -157*P2*z**8/s.Integer(18)+215*z**10/s.Integer(18))
ident(shifted,stated)
print('FIDELITY: Example 6.6 deletion identity verified as an exact polynomial identity')

u,v,w=s.symbols('u v w')
schur = u**3+v**3+w**3-u**2*v-u**2*w-v**2*u-v**2*w-w**2*u-w**2*v+3*u*v*w
ident(schur,(u-v)**2*(u+v-w)+w*(u-w)*(v-w))
print('POSITIVITY: sorted-square Robinson decomposition verified exactly')

for n in range(1,8):
    xs=s.symbols('x:'+str(n))
    h=s.Poly(hp(xs),*xs)
    ident(h.as_expr(),ht(xs))
    assert all(sum(a)==18 and all(e%2==0 for e in a) for a,c in h.terms() if c)
    ident(h.as_expr().subs(dict.fromkeys(xs,1)),0)
    print('IDENTITY: H power sums = Robinson triple sum, n='+str(n)+', terms='+str(sum(1 for a,c in h.terms() if c)))
# Independent rational boundary samples; the decomposition above supplies the general sign argument.
for xyz in product(range(-3,4),repeat=3):
    assert rob(*xyz)>=0
print('SAMPLES: Robinson >= 0 on all 343 integer triples in [-3,3]^3')


def moment(alpha):
    assert sum(alpha)==18
    if any(e%2 for e in alpha):
        return Q(0)
    return MU[tuple(sorted((e//2 for e in alpha),reverse=True))]


def matrix(basis):
    return [[moment(tuple(a+b for a,b in zip(x,y))) for y in basis] for x in basis]


def ldl(M):
    n=len(M)
    L=[[Q(int(i==j)) for j in range(n)] for i in range(n)]
    D=[]
    for j in range(n):
        d=M[j][j]-sum((L[j][r]**2*D[r] for r in range(j)),Q(0))
        assert d>0
        D.append(d)
        for i in range(j+1,n):
            L[i][j]=(M[i][j]-sum((L[i][r]*D[r]*L[j][r] for r in range(j)),Q(0)))/d
    for i in range(n):
        for j in range(n):
            assert M[i][j]==sum((L[i][r]*D[r]*L[j][r] for r in range(n)),Q(0))
    return L,D

B100=[(2*a+1,2*b,2*(4-a-b)) for a in range(5) for b in range(5-a)]
B010=[(b,a,c) for a,b,c in B100]
B001=[(c,b,a) for a,b,c in B100]
B111=[(2*a+1,2*b+1,2*(3-a-b)+1) for a in range(4) for b in range(4-a)]
blocks=[B100,B010,B001,B111]
basis=sum(blocks,[])
assert len(basis)==len(set(basis))==55
assert set(basis)=={(a,b,9-a-b) for a in range(10) for b in range(10-a)}
for i in range(4):
    for j in range(i):
        assert all(moment(tuple(a+b for a,b in zip(x,y)))==0 for x in blocks[i] for y in blocks[j])
assert matrix(B100)==matrix(B010)==matrix(B001)
print('GRAM: complete degree-9 basis, 55 monomials; parity sizes 15,15,15,10; cross blocks zero')
all_entries=[]
for bi,B in enumerate([B100,B111]):
    M=matrix(B)
    L,D=ldl(M)
    floors=[(DEN*d.numerator)//d.denominator for d in D]
    assert floors==EXPECTED_FLOORS[bi]
    all_entries += [v for row in L for v in row]+D
    print('LDL: block='+str(len(B))+', reconstruction exact, pivots all positive')
    print('FLOOR_1E9_PIVOTS: '+json.dumps(floors))
print('COST_DATA: '+json.dumps({
    'representative_matrix_entries':15**2+10**2,
    'representative_LDL_entries':15**2+15+10**2+10,
    'max_numerator_bits':max(abs(q.numerator).bit_length() for q in all_entries),
    'max_denominator_bits':max(q.denominator.bit_length() for q in all_entries),
    'LDL_rational_text_bytes':sum(len(str(q)) for q in all_entries)}))

x,y,z=s.symbols('x y z')
def functional(poly):
    return sum((Q(int(c.p),int(c.q))*moment(a) for a,c in s.Poly(poly,x,y,z).terms()),Q(0))
h3=s.expand(hp((x,y,z)))
assert functional(h3)==Q(-3231,DEN)
assert functional(p((x,y,z),18))==3
E=s.symbols('epsilon')
# The affine expression is exact; its unique zero is the family endpoint.
value=3*E-s.Rational(3231,DEN)
assert s.solve(value,E)==[s.Rational(1077,DEN)]
eps=Q(1,10**7)
f3=h3+s.Rational(eps.numerator,eps.denominator)*p((x,y,z),18)
assert functional(f3)==Q(-2931,DEN)
print('SEPARATION: L(H)=-3231/1000000000; L(p18)=3')
print('FAMILY: L(F_epsilon)=3*epsilon-3231/1000000000 < 0 iff epsilon<1077/1000000000')
print('WITNESS: epsilon=1/10000000; L(F^(3))=-2931/1000000000')

# Exact tail substitution with N=square^2 avoids introducing approximate square roots.
t=s.symbols('t')
for m in (1,2,3):
    head=(x,y,z)[:m]
    for square in (1,2,3):
        N=square**2
        tail=(t/s.Integer(square),)*N
        for j in range(2,19):
            ident(p(head+tail,j),p(head,j)+s.Rational(square)**(2-j)*t**j)
print('TAIL: exact substitution for m=1,2,3; N=1,4,9; powers j=2,...,18')
for k in range(5):
    P=s.Poly((p((x,y,z),2)+t**2)**k*f3,t)
    ident(P.coeff_monomial(t**(2*k)),f3)
    assert P.degree()==2*k
print('EXTRACTION: coefficient of t^(2k) equals F^(3), k=0,...,4')
# On the n-coordinate flat unit vector, H=0 and p18=n/n^9=n^-8.
for n in range(1,33):
    assert Q(n,n**9)==Q(1,n**8)
print('FLAT_SPHERE: exact p18=n^-8 for n=1,...,32, after H vanishes')
print('ALL_CHECKS_OK')
