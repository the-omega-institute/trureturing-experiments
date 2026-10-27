#!/usr/bin/env python3
"""Exact checks for the Section V question in arXiv:2610.01940v1.

Python 3 and SymPy 1.14; all constants are integers, rationals, or exact algebraic
and symbolic expressions. No numerical approximations are used. Finite circuit
checks support, but do not replace, the unrestricted circuit induction.
"""
from itertools import product
import sympy as s

I, Q = s.I, s.Rational
E2, E4 = s.eye(2), s.eye(4)
Had = s.Matrix([[1, 1], [1, -1]]) / s.sqrt(2)
X = s.Matrix([[0, 1], [1, 0]])
phase = s.diag(1, I)


def clean(x):
    return s.simplify(s.expand_complex(x))


def matrix_equal(a, b):
    assert a.shape == b.shape
    assert all(s.simplify(x) == 0 for x in a - b), (a, b)
    assert not a.atoms(s.Float) and not b.atoms(s.Float)


def laurent_support(expr, z):
    expr = s.expand(expr)
    powers = set()
    for term in s.Add.make_args(expr):
        if term == 0:
            continue
        k = term.as_powers_dict().get(z, s.Integer(0))
        assert k.is_Integer
        assert not (term / z ** k).has(z)
        powers.add(int(k))
    return powers


# Fidelity: Appendix C, proof of Lemma `lem:query-polynomial`, states
# A(z) = sum_{k=-N}^N A_k z^k, and its channel entries have support in [-2N,2N].
# Reproduce this result exactly for every sign word through N=3, nontrivial
# oracle-independent mixing gates, and the source's hard family at delta=pi/4.
z = s.Symbol('z', nonzero=True)
a = (1 + I) / s.sqrt(2)
assert 0 < s.pi / 4 < 1
O = {1: s.diag(a, z, 1, 1), -1: s.diag(s.conjugate(a), 1 / z, 1, 1)}
W = [s.kronecker_product(Had, Had),
     s.kronecker_product(phase, Had),
     s.kronecker_product(X, phase)]
for w in W:
    matrix_equal(w.H * w, E4)


def star_laurent(expr):
    return s.conjugate(expr).subs(s.conjugate(z), 1 / z)


fidelity_count = 0
for n in range(4):
    for signs in product((-1, 1), repeat=n):
        A = W[0]
        for j, sign in enumerate(signs):
            A = (W[(j + 1) % len(W)] * O[sign] * A).applyfunc(s.expand)
        for entry in A:
            assert all(-n <= k <= n for k in laurent_support(entry, z))
        for row in range(4):
            for col in range(4):
                entry = s.expand(A[row, col] * star_laurent(A[(row + 1) % 4, (col + 1) % 4]))
                assert all(-2 * n <= k <= 2 * n for k in laurent_support(entry, z))
        fidelity_count += 1
print('FIDELITY: source Laurent amplitude/channel support bounds; all 15 sign words for N=0..3 OK')
assert fidelity_count == 15

# Exact oracle pair, prescribed Hermitian logarithms, and shared spectral promise.
U, V = s.diag(I, -I), s.diag(-I, I)
HU, HV = s.diag(s.pi / 2, 3 * s.pi / 2), s.diag(3 * s.pi / 2, s.pi / 2)
for H, B in [(HU, U), (HV, V)]:
    matrix_equal(H.H, H)
    matrix_equal((I * H).exp().applyfunc(clean), B)
    matrix_equal(B.H * B, E2)
    assert len(H.eigenvals()) == 2
    assert B[0, 0] != B[1, 1]
    for eigenvalue in H.eigenvals():
        assert 1 < eigenvalue < 2 * s.pi
matrix_equal(V, -U)
print('PROMISE: nonscalar U,V; exp(i H)=oracle; both spectra lie in (1,2*pi), hence in [delta,2*pi) for every 0<delta<1')

# Finite independent confirmation of the channel induction, with and without
# a qubit of workspace. The outline supplies the arbitrary-N argument.
phase_cases = 0
for workspace in (1, 2):
    size = 2 * workspace
    gates = ([Had, phase, X] if workspace == 1 else W)
    for n in range(5):
        for signs in product((-1, 1), repeat=n):
            outputs = []
            for B in (U, V):
                A = gates[0]
                for j, sign in enumerate(signs):
                    query = B if sign == 1 else B.H
                    A = gates[(j + 1) % len(gates)] * s.kronecker_product(s.eye(workspace), query) * A
                outputs.append(A.applyfunc(s.simplify))
            AU, AV = outputs
            matrix_equal(AV, (-1) ** n * AU)
            matrix_equal(s.kronecker_product(AU, s.conjugate(AU)),
                         s.kronecker_product(AV, s.conjugate(AV)))
            phase_cases += 1
assert phase_cases == 62
print('INDISTINGUISHABILITY: 62 exact circuits; both signs, N=0..4, workspace dimensions 1 and 2 OK')

# All-real-t symbolic target outputs and separation identity.
t = s.Symbol('t', real=True)
rho = s.ones(2) / 2
Ut = s.diag(s.exp(I * s.pi * t / 2), s.exp(3 * I * s.pi * t / 2))
Vt = s.diag(s.exp(3 * I * s.pi * t / 2), s.exp(I * s.pi * t / 2))
R = (Ut * rho * Ut.H).applyfunc(clean)
S = (Vt * rho * Vt.H).applyfunc(clean)
D = s.Matrix([[0, -I * s.sin(s.pi * t)], [I * s.sin(s.pi * t), 0]])
matrix_equal(R - S, D)
matrix_equal(D.H, D)
matrix_equal(D.H * D, s.sin(s.pi * t) ** 2 * E2)
lam = s.Symbol('lambda')
assert s.simplify(D.charpoly(lam).as_expr() - (lam ** 2 - s.sin(s.pi * t) ** 2)) == 0
print('FAMILY: symbolic D(t)^*D(t)=sin(pi*t)^2 I; singular values |sin(pi*t)|, trace norm 2|sin(pi*t)|')

# Explicit rational-t instances, including the preregistered contradiction.
for value in (Q(1, 6), Q(1, 4), Q(1, 3), Q(1, 2), Q(2, 3), Q(3, 4), Q(5, 6)):
    difference = D.subs(t, value).applyfunc(clean)
    eigenvalues = difference.eigenvals()
    exact_norm = s.simplify(sum(s.Abs(e) * mult for e, mult in eigenvalues.items()))
    assert s.simplify(exact_norm - 2 * s.sin(s.pi * value)) == 0
Rhalf, Shalf = R.subs(t, Q(1, 2)), S.subs(t, Q(1, 2))
for M in (Rhalf, Shalf):
    matrix_equal(M * M, M)
    assert M.trace() == 1
matrix_equal(Rhalf * Shalf, s.zeros(2))
epsilon = Q(1, 4)
assert 2 > 2 * epsilon
print('FIXED WITNESS: delta=1/2, t=1/2, epsilon=1/4; orthogonal rank-one targets, norm separation 2 > 2*epsilon=1/2')
print('BOUNDARY: at t=0 or t=1 the target-channel difference is zero')
for value in (0, 1):
    matrix_equal(D.subs(t, value), s.zeros(2))
print('LIMIT: finite circuit checks and symbolic identities are not a Lean proof; arbitrary-N induction and diamond error inference are in the draft outline')
print('ALL_CHECKS_OK')
