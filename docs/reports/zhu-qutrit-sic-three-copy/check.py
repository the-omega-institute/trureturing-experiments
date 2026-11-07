#!/usr/bin/env python3
"""Exact checks for Zhu, arXiv:2111.02841v2, Conjecture 1.

Python 3.12, SymPy 1.14.0. All equalities are over rational Laurent
polynomials or exact algebraic expressions. No floating-point evaluations.
This checks the proposed family calculation, not the all-SIC classification
or a Lean proof of the conjecture. Source: DOI 10.1103/PRXQuantum.3.030306.
"""
from collections import Counter
from itertools import combinations, permutations, product
import sympy as s

R = s.Rational
t, w, x = s.symbols('t w x')
K = s.QQ.frac_field(t)
modulus = s.Poly(w*w+w+1, w, domain=K)

def red(e):
    return s.Poly(s.expand(e), w, domain=K).rem(modulus).as_expr().expand()

def mat(M):
    return M.applyfunc(red)

def adj(M):
    return mat(M.T.applyfunc(lambda e: e.subs({t:1/t, w:w*w}, simultaneous=True)))

def eq(A, B):
    if isinstance(A, s.MatrixBase):
        assert mat(A-B) == s.zeros(*A.shape), (A, B)
    else:
        assert red(A-B) == 0, (A, B)

def cp(M):
    return s.Poly(M.charpoly(x).as_expr(), x).all_coeffs()

def cp_eq(M, polynomial):
    a, b = cp(M), s.Poly(polynomial, x).all_coeffs()
    assert len(a) == len(b)
    for u,v in zip(a,b): eq(u,v)

I = s.eye(3)
X = s.Matrix([[0,0,1],[1,0,0],[0,1,0]])
Z = s.diag(1,w,w*w)
v = s.Matrix([0,1,-t])
vecs = [mat(X**j*Z**k*v) for j,k in product(range(3),repeat=2)]
P = [mat(a*adj(a)/2) for a in vecs]
for A in P:
    eq(adj(A),A); eq(A*A,A); eq(s.trace(A),1)
eq(sum(P,s.zeros(3)),3*I)
for i,j in combinations(range(9),2): eq(s.trace(P[i]*P[j]),R(1,4))
print('SIC projectors: 9 idempotents, 36 pair overlaps, sum=3I (symbolic t, t*=1/t)')

# Source Eq. (47), with trace-one arguments. Scalar terms multiply I.
def source_Q(A,B,C):
    f12=red(s.trace(A*B));f23=red(s.trace(B*C));f31=red(s.trace(C*A))
    scalar=1+f12+f23+f31+red(s.trace(A*B*C))+red(s.trace(A*C*B))
    return mat(scalar*I+(1+f23)*A+(1+f31)*B+(1+f12)*C+
               A*B+B*A+A*C+C*A+B*C+C*B+
               A*B*C+A*C*B+B*C*A+B*A*C+C*A*B+C*B*A)

# Independent Eq. (36) contraction: sum over S_4, without using Eq. (47).
# Each closed permutation cycle contributes a trace; the cycle through output
# slot 3 contributes the ordered product of its other slots.
def permutation_Q(A,B,C):
    matrices=[A,B,C]; answer=s.zeros(3)
    for p in permutations(range(4)):
        cycles=[]; seen=set()
        for j in range(4):
            if j in seen: continue
            cy=[]; k=j
            while k not in seen:
                seen.add(k);cy.append(k);k=p[k]
            cycles.append(cy)
        scalar=s.Integer(1); opened=I
        for cy in cycles:
            if 3 in cy:
                n=cy.index(3); order=cy[n+1:]+cy[:n]
                for k in order: opened=opened*matrices[k]
            else:
                m=I
                for k in cy:m=m*matrices[k]
                scalar*=s.trace(m)
        answer+=scalar*opened
    return mat(answer)
# Exact generic rational matrices test every term of the contraction formula.
A=s.Matrix([[1,2,0],[0,1,3],[4,0,2]])
B=s.Matrix([[2,0,1],[1,3,0],[0,2,1]])
C=s.Matrix([[1,0,2],[3,2,1],[0,1,4]])
# Generic formula includes arbitrary traces, unlike the trace-one source_Q.
def generic_Q(A,B,C):
    a,b,c=map(s.trace,(A,B,C))
    return mat((a*b*c+s.trace(A*B)*c+s.trace(B*C)*a+s.trace(C*A)*b+
                s.trace(A*B*C)+s.trace(A*C*B))*I+
               b*c*A+c*a*B+a*b*C+s.trace(B*C)*A+s.trace(C*A)*B+s.trace(A*B)*C+
               c*(A*B+B*A)+b*(A*C+C*A)+a*(B*C+C*B)+
               A*B*C+A*C*B+B*C*A+B*A*C+C*A*B+C*B*A)
eq(permutation_Q(A,B,C),generic_Q(A,B,C))
print('Source Eq. (43): independent 24-permutation contraction agrees for exact integer matrices')

# Reproduce source two-copy SIC fidelity: Eqs. (37), (50), (65).
for i,j in combinations(range(9),2):
    A,B=P[i],P[j]
    Q2=mat(R(5,4)*I+A+B+A*B+B*A)
    cp_eq(Q2,(x-R(7,2))*(x-R(5,4))*(x-R(3,2)))
f2=(9*6+72*R(7,2))/(9*3*4*5)
assert f2==R(17,30)
assert f2==(2*3*4+2*(1+(3-1)*s.sqrt(3+1)))/(3*4*5)
print('Source two-copy qutrit SIC fidelity: 17/30 (all 81 ordered outcomes)')

# All 84 distinct unordered triples: matrix reduction and Newton cubic.
classes=Counter()
for i,j,k in combinations(range(9),3):
    A,B,C=P[i],P[j],P[k];S=mat(A+B+C)
    triple=red(s.trace(A*B*C));r=red((triple+adj(s.Matrix([[triple]]))[0,0])/2)
    classes[r]+=1
    Q=source_Q(A,B,C)
    eq(Q,(R(7,4)+2*r)*I+S**3-S**2+R(3,4)*S)
    cp_eq(S,(x-1)**3-R(3,4)*(x-1)-2*r)
    eq(Q,permutation_Q(A,B,C))
expected=Counter({-R(1,8):3,R(1,16):54})
for r in [-(t**3+t**-3)/16,
          ((1+w)*t**3-w*t**-3)/16,
          (-w*t**3+(1+w)*t**-3)/16]:expected[red(r)]+=9
assert classes==expected,(classes,expected)
print('Distinct triples: 84 matrix/cubic/contraction checks; phase multiplicities 3,54,9,9,9')

# Repeated outcomes: exact characteristic polynomials; no eigensolver rounding.
for i in range(9):
    cp_eq(source_Q(P[i],P[i],P[i]),(x-24)*(x-6)**2)
for i,j in permutations(range(9),2):
    cp_eq(source_Q(P[i],P[i],P[j]),(x-3)*(x*x-15*x+R(81,2)))
counts=Counter(len(set(a)) for a in product(range(9),repeat=3))
assert counts=={1:9,2:216,3:504}
print('Ordered outcomes: 9 all equal; 216 two equal; 504 distinct; maximum roots 24, (15+3*sqrt(7))/2')

# Cubic spectral and harmonic identities as rational Laurent identities.
z=s.symbols('z',nonzero=True)
cos=lambda k:(z**k+z**(-k))/2
r=cos(3)/8;lam=1+cos(1)
assert s.expand((lam-1)**3-R(3,4)*(lam-1)-2*r)==0
q=lambda y:R(7,2)+(y*y+y**-2)/2+(y**3+y**-3)/4+R(5,4)*(y+y**-1)
assert s.expand(R(7,4)+2*r+lam**3-lam**2+R(3,4)*lam-q(z))==0
h=x**3-x*x+R(3,4)*x
assert s.expand(s.diff(h,x)-(3*(x-R(1,3))**2+R(5,12)))==0
b=s.symbols('b',nonzero=True)
AA=R(5,2)*(1+b+b**-1);BB=1+b*b+b**-2
harmonic=s.expand(q(b*z)+q(z)+q(b/z)-R(21,2)-AA*(z+1/z)/2-BB*(z*z+z**-2)/2)
num=s.fraction(s.cancel(harmonic))[0]
assert s.rem(num,b**6+b**3+1,b)==0
print('Spectral cubic, increasing polynomial h, and third-harmonic cancellation: exact identities')

# Normalization and closed formula in the real algebraic field c=cos(pi/9).
c,g=s.symbols('c g')
modc=8*c**3-6*c-1
rc=lambda e:s.rem(s.expand(e),modc,c).expand()
cos2=2*c*c-1
Acoef=R(5,2)*(1+2*cos2)
cos4=rc(8*c**4-8*c*c+1)
Bcoef=1+2*cos4
Ccoef=R(91,240)+g/30+cos2/30+c/12
qv=R(15,4)+cos2+R(5,2)*c
constant=(9*24+216*(15+3*g)/2+6*(3*R(15,4)+54*qv+9*R(21,2)))/9720
assert rc(constant-Ccoef)==0
assert R(6*9,9720)==R(1,180)
assert rc(cos4-(1-c)*(2*c+1))==0
assert rc(Acoef-R(5,2)*(4*c*c-1))==0
assert rc(Bcoef-(1+2*(1-c)*(2*c+1)))==0
# c is the largest root, isolated by rational endpoints; no float evidence.
intervals=s.Poly(modc,c).intervals(eps=R(1,1000))
assert len(intervals)==3 and intervals[-1][0][0]>R(1,2) and intervals[-1][0][1]<1
Fzero=Ccoef+(Acoef*c+Bcoef*cos2)/180
Fend=Ccoef+(Acoef+Bcoef)/180
assert rc(Fend-Fzero-(1-c)*(Acoef+2*Bcoef*(1+c))/180)==0
print('Three-copy formula: C+[A*cos(u)+B*cos(2*u)]/180, u=pi/9-phi')
print('C=91/240+sqrt(7)/30+cos(2*pi/9)/30+cos(pi/9)/12')
print('A=(5/2)*(1+2*cos(2*pi/9)); B=1+2*cos(4*pi/9)')
print('cos(pi/9) exact isolating interval:',intervals[-1][0])
print('Endpoint difference=(1-c)*(A+2*B*(1+c))/180; A>0 and B>0 for 1/2<c<1')
print('Scope: exact family identities; classification and analytic branch/order bridges remain proof obligations')
print('ALL_CHECKS_OK')
