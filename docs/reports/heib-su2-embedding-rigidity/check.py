#!/usr/bin/env python3
"""Exact finite checks for Heib's compact-embedding Conjecture 33.

These checks reproduce source identities and test bounded families. They do not
prove the unbounded conjecture, PBW, or the abelian-centralizer theorem.
"""
from functools import lru_cache
from itertools import product
from math import factorial
import sympy as s
from sympy import QQ_I

D = QQ_I
ZERO, ONE = D.zero, D.one
I = D.convert(s.I)


def clean(p):
    return {w: c for w, c in p.items() if c != ZERO}


def add(*ps):
    ans = {}
    for p in ps:
        for w, c in p.items():
            ans[w] = ans.get(w, ZERO) + c
    return clean(ans)


def scale(c, p):
    c = D.convert(c)
    return clean({w: c * a for w, a in p.items()})


@lru_cache(None)
def normal(word):
    """Ordered PBW words X^a Y^b Z^c; apply adjacent inversion rules."""
    for k in range(len(word) - 1):
        a, b = word[k:k + 2]
        if a > b:
            swapped = word[:k] + (b, a) + word[k + 2:]
            # YX=XY-Z; ZX=XZ+Y; ZY=YZ-X.
            gen, sign = {(1, 0): (2, -1), (2, 0): (1, 1),
                         (2, 1): (0, -1)}[(a, b)]
            lower = word[:k] + (gen,) + word[k + 2:]
            ans = dict(normal(swapped))
            for w, c in normal(lower):
                ans[w] = ans.get(w, ZERO) + sign * c
            return tuple(sorted(clean(ans).items()))
    return ((word, ONE),)


def mul(p, q):
    ans = {}
    for u, a in p.items():
        for v, b in q.items():
            for w, c in normal(u + v):
                ans[w] = ans.get(w, ZERO) + a * b * c
    return clean(ans)


def power(p, n):
    ans = {(): ONE}
    for _ in range(n):
        ans = mul(ans, p)
    return ans


def comm(p, q):
    return add(mul(p, q), scale(-1, mul(q, p)))


def dagger(p):
    ans = {}
    for word, coeff in p.items():
        coeff = D.convert(s.conjugate(D.to_sympy(coeff))) * (-1)**len(word)
        for w, c in normal(tuple(reversed(word))):
            ans[w] = ans.get(w, ZERO) + coeff * c
    return clean(ans)


X, Y, Z = ({(j,): ONE} for j in range(3))
Q = scale(-1, add(*(mul(a, a) for a in (X, Y, Z))))
assert comm(X, Y) == Z and comm(Y, Z) == X and comm(Z, X) == Y
assert all(not comm(a, Q) for a in (X, Y, Z))
assert dagger(Q) == Q
# Reproduce Example 26 in the full enveloping algebra, not a central quotient.
E = add(X, scale(-s.I, Y))
F = add(scale(-1, X), scale(-s.I, Y))
H = scale(-2*s.I, Z)
Ep = add(E, scale(s.Rational(1, 4), add(mul(F, H), mul(H, F))),
         scale(-s.Rational(1, 4), power(F, 3)))
Fp, Hp = F, add(H, scale(-1, power(F, 2)))
assert comm(Hp, Ep) == scale(2, Ep)
assert comm(Hp, Fp) == scale(-2, Fp)
assert comm(Ep, Fp) == Hp
Bp = (scale(s.Rational(1, 2), add(Ep, scale(-1, Fp))),
      scale(s.I/2, add(Ep, Fp)), scale(s.I/2, Hp))
assert comm(Bp[0], Bp[1]) == Bp[2]
assert comm(Bp[1], Bp[2]) == Bp[0]
assert comm(Bp[2], Bp[0]) == Bp[1]
assert any(dagger(b) != scale(-1, b) for b in Bp)
assert scale(-1, add(*(mul(b, b) for b in Bp))) == Q
print('FIDELITY: Example 26 nonlinear sl2 brackets hold in exact PBW arithmetic;')
print('          its compact triple fails the fixed dagger, while its Casimir equals Q.')

# Test reduction confluence on every word of length at most six by selecting each
# possible inversion as the first reduction, then using the normal reducer.
words_tested = 0
for length in range(7):
    for w in product(range(3), repeat=length):
        expected = dict(normal(w))
        for k in range(length - 1):
            a, b = w[k:k + 2]
            if a > b:
                gen, sign = {(1, 0): (2, -1), (2, 0): (1, 1),
                             (2, 1): (0, -1)}[(a, b)]
                got = add(dict(normal(w[:k] + (b, a) + w[k + 2:])),
                          scale(sign, dict(normal(w[:k] + (gen,) + w[k + 2:]))))
                assert got == expected
        p = expected
        assert dagger(dagger(p)) == p
        words_tested += 1
print(f'PBW CHECK: {words_tested} words of length 0..6; every first inversion agrees; dagger squared is identity.')


def matrix_zero(a):
    return all(s.expand(x) == 0 for x in a)


def spin(n):
    plus, minus = s.zeros(n + 1), s.zeros(n + 1)
    for k in range(n + 1):
        if k < n:
            plus[k + 1, k] = n - k
        if k > 0:
            minus[k - 1, k] = k
    return (-s.I*s.Rational(1, 2)*(plus + minus),
            s.Rational(1, 2)*(plus - minus),
            s.diag(*(s.I*(k - s.Rational(n, 2)) for k in range(n + 1))))


def evaluate(p, mats):
    ans = s.zeros(mats[0].rows)
    for word, coeff in p.items():
        a = s.eye(mats[0].rows)
        for j in word:
            a = a*mats[j]
        ans += D.to_sympy(coeff)*a
    return ans.applyfunc(s.expand)


for n in range(9):
    mats = spin(n)
    gram = s.diag(*(factorial(k)*factorial(n-k) for k in range(n+1)))
    for a in mats:
        assert matrix_zero(a.conjugate().T*gram + gram*a)
    for i in range(3):
        assert matrix_zero(mats[i]*mats[(i+1)%3] - mats[(i+1)%3]*mats[i] - mats[(i+2)%3])
    cn = s.Rational(n*(n+2), 4)
    assert matrix_zero(evaluate(Q, mats) - cn*s.eye(n+1))
    assert matrix_zero(evaluate(scale(-1, add(*(mul(b, b) for b in Bp))), mats)
                       - cn*s.eye(n+1))
    # Faithful PBW computation is cross-checked against exact matrices.
    for a, b in product((X, Y, Z, Q), repeat=2):
        assert matrix_zero(evaluate(mul(a, b), mats) - evaluate(a, mats)*evaluate(b, mats))
print('FIDELITY: Lemma 31 polynomial spin matrices, n=0..8, satisfy cyclic brackets')
print('          and skew-adjointness for the source factorial Gram matrix.')
print('SETTLEMENT CHECK: their positive Casimir is n(n+2)/4, n=0..8.')


def partitions(total, least=1):
    if total == 0:
        yield ()
    for first in range(least, total+1):
        for tail in partitions(total-first, first):
            yield (first,)+tail


parts_checked = 0
for n in range(9):
    cn = s.Rational(n*(n+2), 4)
    for dims in partitions(n+1):
        values = [s.Rational(d*d-1, 4) for d in dims]
        assert all(0 <= v <= cn for v in values)
        if len(set(values)) == 1:
            assert sum(dims) == n+1 and 0 <= values[0] <= cn
        parts_checked += 1
print(f'SETTLEMENT CHECK: all {parts_checked} irreducible-dimension partitions in dimensions 1..9 obey the Casimir bound.')


def cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def det(rows):
    c = cross(rows[1], rows[2])
    return sum(rows[0][i]*c[i] for i in range(3))


valid, injective = 0, []
for flat in product((-1, 0, 1), repeat=9):
    rows = tuple(flat[3*j:3*j+3] for j in range(3))
    if all(cross(rows[j], rows[(j+1)%3]) == rows[(j+2)%3] for j in range(3)):
        # All purely imaginary scalar offsets in {-i,0,i}^3 are tested too.
        for offsets in product((-1, 0, 1), repeat=3):
            if all(offsets[(j+2)%3] == 0 for j in range(3)):
                valid += 1
                if det(rows) != 0:
                    assert det(rows) == 1
                    assert all(sum(rows[i][k]*rows[j][k] for k in range(3)) == int(i == j)
                               for i, j in product(range(3), repeat=2))
                    injective.append(rows)
                else:
                    assert flat == (0,)*9
assert valid == 25 and len(injective) == 24
print('SETTLEMENT CHECK: all 531441 affine triples with real generator coefficients')
print('          in {-1,0,1} and imaginary scalar coefficients in {-i,0,i}:')
print('          exactly 25 Lie maps, 24 embeddings; every embedding has the canonical image.')

# Radial nonlinear family B_i=f(Q)X_i. Q is central above, so the bracket
# condition is f^2=f, an exact polynomial identity (not evaluation at spins).
t = s.symbols('t')
radial_valid = []
for coeffs in product(range(-2, 3), repeat=3):
    f = s.Poly(sum(c*t**k for k, c in enumerate(coeffs)), t, domain=s.QQ)
    fQ = add({(): D.convert(coeffs[0])}, scale(coeffs[1], Q),
             scale(coeffs[2], power(Q, 2)))
    defect = add(mul(fQ, fQ), scale(-1, fQ))
    actual_brackets_hold = all(not mul(defect, a) for a in (X, Y, Z))
    assert actual_brackets_hold == (f*f == f)
    if actual_brackets_hold:
        radial_valid.append(coeffs)
assert radial_valid == [(0,0,0), (1,0,0)]
print('SETTLEMENT CHECK: all 125 central-polynomial radial triples of polynomial degree <=2,')
print('          coefficients -2..2: only zero and the canonical embedding satisfy the brackets.')
# Negative controls protect the bracket and normalization tests.
assert comm(scale(2, X), Y) != Z
assert not matrix_zero(evaluate(Q, spin(1)) - s.Rational(3, 2)*s.eye(2))
print('LIMIT: finite checks support the route; the full unbounded conjecture is not proved by this program.')
print('ALL_CHECKS_OK')
