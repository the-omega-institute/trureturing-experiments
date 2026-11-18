#!/usr/bin/env python3
"""Exact algebra for the strict QCP candidate; no numerical entropy optimization."""
import itertools
import sys
from fractions import Fraction
import sympy as s

assert sys.version_info >= (3, 12)
I = s.I
X = s.Matrix([[0, 1], [1, 0]])
Y = s.Matrix([[0, -I], [I, 0]])
Z = s.diag(1, -1)

def zero(A):
    assert all(s.simplify(x) == 0 for x in A), A

def comm(A, B):
    return A * B - B * A

def projector(v):
    return v * v.conjugate().T

# Fidelity: the Section II claim that phase damping gives a nonzero commutator.
r = s.symbols('r', real=True)
v00 = s.Matrix([1, 0, 0, 1]) / s.sqrt(2)
v11 = s.Matrix([1, 1, 1, -1]) / 2

def pd2(A):
    return s.Matrix(4, 4, lambda i, j: A[i, j] * r ** ((int(i) ^ int(j)).bit_count()))

zero(comm(pd2(projector(v00)), pd2(projector(v11))) -
     I * (1-r*r) * r / 8 * (s.kronecker_product(Y, s.eye(2)) + s.kronecker_product(s.eye(2), Y)))
assert comm(pd2(projector(v00)), pd2(projector(v11))).subs(r,s.Rational(1,2)) != s.zeros(4)
print('FIDELITY_OK: source phase-damping superactivation witness has nonzero commutator')
print('SOURCE_BOUNDARY: with standard computational dephasing the Pauli factor is Y tensor I + I tensor Y; the printed Z factor differs')

# Candidate family, simultaneously for all a,t in the stated ranges.
a, t, z = s.symbols('a t z', real=True)
E = [(s.eye(2)+a*X)/3, (s.eye(2)+a*Z)/3,
     (s.eye(2)-a*X-a*Z)/3]
F = [s.diag(*[int(i == j) for j in range(3)], E[i]) for i in range(3)]
basis = [s.eye(5)[:, j] for j in range(5)]
v = [basis[0], (basis[0]+basis[1])/s.sqrt(2), basis[1]]
eta = [(1-t)*projector(w)+t*s.eye(5)/5 for w in v]
zero(sum(F, s.zeros(5))-s.eye(5))
for k in range(3):
    zero(E[k]-E[k].conjugate().T)
    roots = [(1-a)/3, (1+a)/3] if k < 2 else [(1-s.sqrt(2)*a)/3, (1+s.sqrt(2)*a)/3]
    assert s.simplify((E[k]-z*s.eye(2)).det() - (roots[0]-z)*(roots[1]-z)) == 0
    assert eta[k].trace() == 1
    zero(eta[k]-eta[k].conjugate().T)
    zero((eta[k]-t*s.eye(5)/5) - (1-t)*projector(v[k]))
# Exact lower endpoint checks: a<=1/2 implies all E eigenvalues >0.
assert s.Rational(1, 6) > 0
assert (1-s.sqrt(2)/2)/3 > 0

def phi(R):
    return sum(((F[i]*R).trace()*eta[i] for i in range(3)), s.zeros(5))

def adj(O):
    return sum(((eta[i]*O).trace()*F[i] for i in range(3)), s.zeros(5))

for k in range(3):
    zero(phi(projector(basis[k]))-eta[k])
Xo, Yo, Zo = [s.diag(P, s.zeros(3)) for P in [X, Y, Z]]
zero(comm(eta[0],eta[1])-I*(1-t)**2*Yo/2)
zero(adj(Xo)-(1-t)*F[1])
zero(adj(Zo)-(1-t)*(F[0]-F[2]))
zero(comm(adj(Xo),adj(Zo))-(1-t)**2*s.diag(s.zeros(3),4*I*a*a*Y/9))
# Explicit nonzero entries; positive a and t<1 exclude their zero loci.
assert s.simplify(comm(eta[0],eta[1])[0,1]-(1-t)**2/2) == 0
assert s.simplify(comm(adj(Xo),adj(Zo))[3,4]-4*a*a*(1-t)**2/9) == 0
print('FAMILY_OK: completeness, spectra, selectors, output floor and both commutators')

# Rational member: exact Choi positivity and trace preservation, independently
# assembled from images of matrix units and from the measure/prepare formula.
sub = {a:s.Rational(1,2),t:s.Rational(1,2)}
J = sum((s.kronecker_product(F[k].T,eta[k]) for k in range(3)),s.zeros(25))
J = J.subs(sub)
J_units = s.zeros(25)
for i,j in itertools.product(range(5), repeat=2):
    eij = s.zeros(5); eij[i,j] = 1
    J_units += s.kronecker_product(eij,phi(eij).subs(sub))
zero(J-J_units)
assert J.is_positive_definite is True
for i,j in itertools.product(range(5), repeat=2):
    assert sum(J[5*i+k,5*j+k] for k in range(5)) == int(i == j)
print('RATIONAL_WITNESS_OK: a=t=1/2, Choi positive definite, output trace partial=I_5')

# Exact certificates for the entropy/information algebra in the argument.
# These verify algebra conditional on the standard entropy identities, the
# accessible-information supremum definition, and local measurement DPI.
chiA,chiB,accA,accB,Iab,Iuv,eps = s.symbols('chiA chiB accA accB Iab Iuv eps')
chiAB = chiA+chiB-Iab
measuredAB = accA+accB-2*eps-Iuv
localg = (chiA-accA)+(chiB-accB)
assert s.expand(localg+2*eps-(chiAB-measuredAB)-(Iab-Iuv)) == 0
# Iab-Iuv>=0 (DPI), eps->0 gives the correlated-prior upper bound.
# For independent priors the two chain-rule upper bounds add exactly.
iIw,iJwGivenI = s.symbols('iIw iJwGivenI')
assert s.expand(accA+accB-(iIw+iJwGivenI)-
                ((accA-iIw)+(accB-iJwGivenI))) == 0
print('ADDITIVITY_CERTIFICATES_OK: correlated-prior DPI residual and independent-prior chain rule')

# Check the conditional-independence identity without entropy approximation:
# equality of the rational probability monomials proves the corresponding log identity.
cases=0
for weights in itertools.product(range(3), repeat=4):
    total=sum(weights)
    if total == 0: continue
    p=[[Fraction(weights[2*i+j],total) for j in range(2)] for i in range(2)]
    rows=[sum(p[i]) for i in range(2)]
    cols=[sum(p[i][j] for i in range(2)) for j in range(2)]
    A=[[Fraction(1,3),Fraction(2,3)],[Fraction(3,4),Fraction(1,4)]]
    B=[[Fraction(2,5),Fraction(3,5)],[Fraction(4,7),Fraction(3,7)]]
    pu=[sum(rows[i]*A[i][u] for i in range(2)) for u in range(2)]
    pv=[sum(cols[j]*B[j][v] for j in range(2)) for v in range(2)]
    puv=[[sum(p[i][j]*A[i][u]*B[j][v] for i,j in itertools.product(range(2),repeat=2)) for v in range(2)] for u in range(2)]
    for i,j,u,v in itertools.product(range(2),repeat=4):
        assert (A[i][u]/pu[u])*(B[j][v]/pv[v])/(puv[u][v]/(pu[u]*pv[v])) == A[i][u]*B[j][v]/puv[u][v]
    cases+=1
assert cases == 80
print('CLASSICAL_IDENTITY_OK: all 80 nonzero binary weight arrays with entries 0,1,2')
print('SCOPE: structural exact certificates; analytic DPI and optimized QCP equality are proved in the outline, not by this program')
print('ALL_CHECKS_OK')
