#!/usr/bin/env python3
"""Exact checks for Aloy et al., arXiv:2312.09278v2, spin-one rotation boxes.

Uses only integer/rational arithmetic and SymPy exact symbolic algebra.
Bounded tests enumerate multisets of rank-one Gaussian-integer Gram matrices;
they do not enumerate all nonnegative real trigonometric polynomials.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement, product
import sympy as s

I = s.I
z = s.Symbol('z', nonzero=True)
phi = s.Matrix([1/z, 1, z])
phi_star = s.Matrix([[z, 1, 1/z]])
J = s.zeros(3)
for i in range(3):
    J[i, 2-i] = 1


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def matrix_zero(mat):
    for x in mat:
        zero(x)


def gram_polynomial(mat):
    return s.expand((phi_star * mat * phi)[0])


def principal_minors(mat):
    return {(ix): s.simplify(mat.extract(ix, ix).det())
            for k in range(1, 4) for ix in combinations(range(3), k)}


# Fidelity: Appendix C.3, proof of Theorem 6, item 2.
# Source order is (1,0,-1); reverse the order to (-1,0,1).
psi = s.Matrix([1, s.sqrt(2), 1])/2
eta = s.Matrix([-1, s.sqrt(2), -1])/2
rho = psi * psi.H
E = eta * eta.H
U = s.diag(1/z, 1, z)
Ustar = s.diag(z, 1, 1/z)
fidelity = s.trace(U * rho * Ustar * E)
zero(fidelity - (2-z-1/z)**2/16)
zero(fidelity - ((1-(z+1/z)/2)**2/4))
zero((psi.H*psi)[0]-1)
matrix_zero(E*E-E)
for effect in (E, s.eye(3)-E):
    assert all(x >= 0 for x in principal_minors(effect).values())
print('FIDELITY_OK: source Theorem 6 appendix realization gives sin(theta/2)^4')

# Symbolic identity on a completely general complex Hermitian 3-by-3 matrix.
d = s.symbols('d0:3', real=True)
r, t, u, v, w, x = s.symbols('r t u v w x', real=True)
Q = s.Matrix([[d[0], r+I*t, u+I*v],
              [r-I*t, d[1], w+I*x],
              [u-I*v, w-I*x, d[2]]])
B = (Q + J*Q.T*J)/2
zero(gram_polynomial(B)-gram_polynomial(Q))
zero(B[0,1]-B[1,2])
zero(B[0,0]-B[2,2])
zero(B[0,2]-Q[0,2])
# Normalization: trace=1, second Fourier coefficient=0, first coefficient=0.
C = B.subs({d[1]: 1-d[0]-d[2], u: 0, v: 0, w: -r, x: -t})
matrix_zero(C-s.diag((d[0]+d[2])/2, 1-d[0]-d[2], (d[0]+d[2])/2))
zero(gram_polynomial(C)-1)
print('SYMBOLIC_OK: reversal-transpose preserves the polynomial and diagonalizes every normalized sum')

# All seven nonempty supports, with arbitrary positive entries on that support.
cs = s.symbols('c0:3', positive=True)
m = s.Symbol('m', positive=True, integer=True)
for bits in product((0, 1), repeat=3):
    if not any(bits):
        continue
    P = s.diag(*bits)
    C = s.diag(*(cs[i] if bits[i] else 0 for i in range(3)))
    D = s.diag(*(s.sqrt(cs[i]) if bits[i] else 0 for i in range(3)))
    M = s.diag(*(1/s.sqrt(cs[i]) if bits[i] else 0 for i in range(3)))
    supported_B = P*B*P
    effect = M*supported_B*M + (s.eye(3)-P)/m
    matrix_zero(M*C*M+s.eye(3)-P-s.eye(3))
    matrix_zero(D*effect*D-supported_B)
    zero(s.trace(U*(D*s.ones(3,1))*(s.ones(1,3)*D)*Ustar*effect)
         -gram_polynomial(supported_B))
print('SUPPORT_OK: all 7 nonempty supports; exact POVM-sum and Born identities')

# Exhaust every unordered family of 3 or 4 factors in {0,+/-1,+/-i}^3,
# modulo a factor's irrelevant global unit phase. Repeated/zero outcomes allowed.
alphabet = (0, 1, -1, I, -I)
vectors = set()
for vals in product(alphabet, repeat=3):
    first = next((a for a in vals if a != 0), None)
    canonical = tuple(s.expand(a*s.conjugate(first)) for a in vals) if first else vals
    vectors.add(canonical)
vectors = sorted(vectors, key=str)
assert len(vectors) == 32
catalog = []
for vals in vectors:
    q = s.Matrix(vals)*s.Matrix(vals).H
    b2 = q+J*q.T*J
    assert all(a >= 0 for a in principal_minors(q).values())
    minors = principal_minors(b2)
    assert all(a >= 0 for a in minors.values())
    zero(gram_polynomial(b2)/2-gram_polynomial(q))
    sig = tuple(int(a) for y in (q[0,1]+q[1,2], q[0,2])
                for a in (s.re(y), s.im(y)))
    flat = tuple((int(s.re(y)),int(s.im(y))) for y in b2)
    catalog.append((int(s.trace(q)), sig, flat, minors))

counts, ranks = Counter(), Counter()
visited = Counter()
for outcomes in (3, 4):
    for family in combinations_with_replacement(range(32), outcomes):
        visited[outcomes] += 1
        total = sum(catalog[a][0] for a in family)
        if not total or any(sum(catalog[a][1][j] for a in family) for j in range(4)):
            continue
        summed = tuple(tuple(sum(catalog[a][2][j][k] for a in family)
                             for k in range(2)) for j in range(9))
        assert all(summed[3*i+j] == (0,0) for i in range(3) for j in range(3) if i != j)
        c = tuple(Fraction(summed[4*i][0], 2*total) for i in range(3))
        assert all(summed[4*i][1] == 0 and c[i] >= 0 for i in range(3))
        assert sum(c) == 1 and c[0] == c[2]
        support = tuple(i for i in range(3) if c[i])
        for a in family:
            _, _, entries, minors = catalog[a]
            # Vanishing rows/columns verify support, including singular sums.
            assert all(entries[3*i+j] == (0,0) for i in range(3) for j in range(3)
                       if i not in support or j not in support)
            # Exact principal-minor certificate for the reconstructed effect.
            # On the support: E=M B M. On its complement: E=I/outcomes.
            for ix in minors:
                active = tuple(i for i in ix if i in support)
                det = Fraction(int(minors[active]), (2*total)**len(active)) if active else Fraction(1)
                for i in active:
                    det /= c[i]
                det /= outcomes**(len(ix)-len(active))
                assert det >= 0
        # Sum E=I follows entry by entry from the diagonal sum and whitening.
        for i in range(3):
            for j in range(3):
                if i == j:
                    value = Fraction(summed[4*i][0],2*total)/c[i] if c[i] else Fraction(1)
                    assert value == 1
                else:
                    assert summed[3*i+j] == (0,0)
        # D E_a D = B_a on support, hence the entire Laurent polynomial agrees;
        # the arbitrary-angle identity, including trace orientation, was checked above.
        counts[outcomes] += 1
        ranks[(outcomes,len(support))] += 1
assert all(counts[n] > 0 for n in (3,4))
assert {rank for n,rank in ranks} == {1,2,3}
for n in (3,4):
    print(f'EXHAUSTIVE_OK: outcomes={n}; all_multisets={visited[n]}; normalized_families={counts[n]}; '
          f'support_ranks='+str({k:ranks[(n,k)] for k in (1,2,3)}))
print('SCOPE: symbolic identities plus the specified finite Gram catalog; no Lean proof or universal factorization certificate')
print('ALL_CHECKS_OK')
