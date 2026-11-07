#!/usr/bin/env python3
"""Exact checks for the finite-field question in arXiv:2611.06305v1, Section 6.

No third-party packages. Weyl products use integral PBW coefficients, checked
against an independent adjacent-word reduction. These finite checks are evidence
for a preregistration, not a proof of its universally quantified conclusion.
"""
from itertools import product, permutations
from math import comb, factorial


class Integers:
    zero, one = 0, 1
    add = staticmethod(lambda a, b: a + b)
    neg = staticmethod(lambda a: -a)
    mul = staticmethod(lambda a, b: a * b)
    cast = staticmethod(int)


class FiniteField:
    """F_p[t]/(monic irreducible modulus), with base-p integer labels."""
    zero, one = 0, 1

    def __init__(self, p, modulus):
        self.p, self.modulus = p, tuple(modulus)
        self.n = len(modulus) - 1
        self.q = p ** self.n
        assert modulus[-1] == 1
        self.elements = range(self.q)
        self.add_table = [[self._add(a, b) for b in self.elements]
                          for a in self.elements]
        self.mul_table = [[self._mul(a, b) for b in self.elements]
                          for a in self.elements]

    def digits(self, a):
        out = []
        for _ in range(self.n):
            a, digit = divmod(a, self.p)
            out.append(digit)
        return out

    def encode(self, digits):
        return sum((a % self.p) * self.p ** i for i, a in enumerate(digits))

    def _add(self, a, b):
        return self.encode([x + y for x, y in zip(self.digits(a), self.digits(b))])

    def _mul(self, a, b):
        c = [0] * (2 * self.n - 1)
        for i, x in enumerate(self.digits(a)):
            for j, y in enumerate(self.digits(b)):
                c[i + j] += x * y
        for k in range(len(c) - 1, self.n - 1, -1):
            t = c[k] % self.p
            for j, coefficient in enumerate(self.modulus):
                c[k - self.n + j] -= t * coefficient
        return self.encode(c[:self.n])

    def add(self, a, b):
        return self.add_table[a][b]

    def neg(self, a):
        return self.encode([-x for x in self.digits(a)])

    def mul(self, a, b):
        return self.mul_table[a][b]

    def cast(self, a):
        return a % self.p

    def power(self, a, n):
        r = 1
        for _ in range(n):
            r = self.mul(r, a)
        return r

    def check_field(self):
        for a in self.elements:
            assert self.add(a, self.neg(a)) == 0
            assert self.power(a, self.q) == a
            if a:
                assert any(self.mul(a, b) == 1 for b in self.elements)
        for a, b, c in product(self.elements, repeat=3):
            assert self.mul(self.mul(a, b), c) == self.mul(a, self.mul(b, c))
            assert self.mul(a, self.add(b, c)) == self.add(self.mul(a, b), self.mul(a, c))


def add(ring, *polys):
    out = {}
    for poly in polys:
        for word, c in poly.items():
            out[word] = ring.add(out.get(word, 0), c)
            if out[word] == 0:
                del out[word]
    return out


def scale(ring, c, poly):
    return {w: ring.mul(c, v) for w, v in poly.items() if ring.mul(c, v) != 0}


def sub(ring, a, b):
    return add(ring, a, scale(ring, ring.neg(1), b))


def weyl_mul(ring, a, b):
    """P^i Q^j * P^k Q^l using QP = PQ - 1."""
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            for t in range(min(j, k) + 1):
                coeff = (-1) ** t * comb(j, t) * factorial(k) // factorial(k - t)
                c = ring.mul(ring.mul(x, y), ring.cast(coeff))
                out = add(ring, out, {(i + k - t, j + l - t): c})
    return out


def free_mul(ring, a, b):
    out = {}
    for u, x in a.items():
        for v, y in b.items():
            out = add(ring, out, {u + v: ring.mul(x, y)})
    return out


def power(ring, mul, a, n, unit):
    out = unit
    for _ in range(n):
        out = mul(ring, out, a)
    return out


def comm(ring, mul, a, b):
    return sub(ring, mul(ring, a, b), mul(ring, b, a))


def gamma3(ring, a, b, c):
    return comm(ring, weyl_mul, comm(ring, weyl_mul, a, b), c)


def standard3(ring, a, b, c):
    out, xs = {}, (a, b, c)
    for perm in permutations(range(3)):
        inversions = sum(perm[i] > perm[j] for i in range(3) for j in range(i + 1, 3))
        term = weyl_mul(ring, weyl_mul(ring, xs[perm[0]], xs[perm[1]]), xs[perm[2]])
        out = add(ring, out, scale(ring, ring.cast((-1) ** inversions), term))
    return out


def t4(ring, a, b, c, d):
    bracket = lambda x, y: comm(ring, weyl_mul, x, y)
    return add(ring,
               weyl_mul(ring, bracket(a, b), bracket(c, d)),
               scale(ring, ring.neg(1), weyl_mul(ring, bracket(a, c), bracket(b, d))),
               weyl_mul(ring, bracket(a, d), bracket(b, c)))


P, Q, ONE = {(1, 0): 1}, {(0, 1): 1}, {(0, 0): 1}


def reduce_word(word):
    """Independent integral reduction of the first QP to PQ - 1."""
    for i in range(len(word) - 1):
        if word[i:i + 2] == (1, 0):
            swapped = word[:i] + (0, 1) + word[i + 2:]
            removed = word[:i] + word[i + 2:]
            return sub(Integers, reduce_word(swapped), reduce_word(removed))
    return {(word.count(0), word.count(1)): 1}


def check_pbw_and_fidelity():
    assert sub(Integers, weyl_mul(Integers, P, Q), weyl_mul(Integers, Q, P)) == ONE
    words = 0
    for n in range(8):
        for word in product(range(2), repeat=n):
            out = ONE
            for x in word:
                out = weyl_mul(Integers, out, (P, Q)[x])
            assert out == reduce_word(word)
            words += 1
    for a, b, c in product((P, Q), repeat=3):
        assert gamma3(Integers, a, b, c) == {}
        assert standard3(Integers, a, b, c) == {}
    for a, b, c, d in product((P, Q), repeat=4):
        assert t4(Integers, a, b, c, d) == {}
    print(f'PBW: independent word reductions agree on {words} words of length <= 7')
    print('FIDELITY: Lemma 6, equations (8)-(9), and (10): all integral basis evaluations vanish (8, 8, 16); multilinearity gives the identities on span(P,Q)')


def univariate_remainder(field, coefficients):
    c = list(coefficients)
    for i in range(len(c) - 1, field.q - 1, -1):
        # Subtract c[i] * t^(i-q) * (t^q - t).
        a, c[i] = c[i], 0
        c[i - field.q + 1] = field.add(c[i - field.q + 1], a)
    return c[:field.q]


def univariate_eval(field, coefficients, x):
    out = 0
    for c in reversed(coefficients):
        out = field.add(field.mul(out, x), c)
    return out


def check_finite_field(field):
    field.check_field()
    vectors = [add(field, scale(field, a, P), scale(field, b, Q))
               for a, b in product(field.elements, repeat=2)]
    assert len({tuple(sorted(v.items())) for v in vectors}) == field.q ** 2
    values = set()
    for a, b, d, e in product(field.elements, repeat=4):
        u = add(field, scale(field, a, P), scale(field, b, Q))
        v = add(field, scale(field, d, P), scale(field, e, Q))
        c = comm(field, weyl_mul, u, v)
        determinant = field.add(field.mul(a, e), field.neg(field.mul(b, d)))
        assert c == ({(0, 0): determinant} if determinant else {})
        values.add(determinant)
        assert sub(field, power(field, weyl_mul, c, field.q, ONE), c) == {}
    assert values == set(field.elements)
    x, y = {(0,): 1}, {(1,): 1}
    c = comm(field, free_mul, x, y)
    h = sub(field, power(field, free_mul, c, field.q, {(): 1}), c)
    degree_two = {w: v for w, v in h.items() if len(w) == 2}
    assert degree_two == scale(field, field.neg(1), c) and degree_two
    assert {len(w) for w in h} == {2, 2 * field.q}
    assert comm(field, weyl_mul, P, Q) == ONE and ONE != {}
    if field.q in (2, 3):
        for a, b, c in product(vectors, repeat=3):
            assert gamma3(field, a, b, c) == {}
            assert standard3(field, a, b, c) == {}
        for a, b, c, d in product(vectors, repeat=4):
            assert t4(field, a, b, c, d) == {}
        print(f'FIDELITY F_{field.q}: exhaustive triples={len(vectors)**3}, quadruples={len(vectors)**4}; generators vanish')
    print(f'FAMILY F_{field.q}: all {field.q**4} pairs in V^2; h_q vanishes; free degree-2 component = -[x,y] != 0; [P,Q]=1 != 0')
    degree = {2: 5, 3: 5, 4: 4, 5: 4}.get(field.q)
    if degree is not None:
        count = 0
        for coefficients in product(field.elements, repeat=degree + 1):
            vanishes = all(univariate_eval(field, coefficients, x) == 0 for x in field.elements)
            divisible = not any(univariate_remainder(field, coefficients))
            assert vanishes == divisible
            count += 1
        print(f'ENVELOPE F_{field.q}: {count} scalar polynomials of degree <= {degree}; vanishing iff divisible by t^q-t; all scalars are commutator values')


def main():
    check_pbw_and_fidelity()
    # Prime fields, then explicit extension fields with irreducible moduli.
    for p, modulus in [(2, [0, 1]), (3, [0, 1]), (2, [1, 1, 1]),
                       (5, [0, 1]), (2, [1, 1, 0, 1]), (3, [1, 0, 1])]:
        check_finite_field(FiniteField(p, modulus))
    print('SCOPE: nonmembership uses the general degree >= 3 filtration argument; finite checks do not certify the infinite ideal or a complete finite-field basis')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
