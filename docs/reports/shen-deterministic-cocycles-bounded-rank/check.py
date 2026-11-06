#!/usr/bin/env python3
"""Exact finite checks for Shen, arXiv:2609.26914v1, Sections 3.3--3.5.

Only integers and Fraction are used. These checks do not prove an infinite
family or implement Shoup's efficient irreducible-polynomial algorithm.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product


def rank_mod(rows, p):
    a = [list(row) for row in rows]
    r = 0
    for j in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][j] % p), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][j], -1, p)
        a[r] = [v * inv % p for v in a[r]]
        for i in range(len(a)):
            if i != r:
                c = a[i][j]
                a[i] = [(v - c * w) % p for v, w in zip(a[i], a[r])]
        r += 1
    return r


def sign(perm):
    return -1 if sum(perm[i] > perm[j] for i in range(3) for j in range(i + 1, 3)) % 2 else 1


PERMS = [(s, sign(s)) for s in permutations(range(3))]


def tensor(p, n, coefficients):
    t = [[[0] * n for _ in range(n)] for _ in range(n)]
    for triple, value in zip(combinations(range(n), 3), coefficients):
        for perm, sg in PERMS:
            i, j, k = [triple[u] for u in perm]
            t[i][j][k] = sg * value % p
    return t


def fidelity():
    expected = {
        (2, 3): {22: 1, 64: 1},
        (2, 4): {88: 15, 256: 1},
        (2, 5): {184: 868, 352: 155, 1024: 1},
        (3, 3): {105: 2, 729: 1},
        (3, 4): {945: 80, 6561: 1},
        (5, 3): {745: 4, 15625: 1},
    }
    for (p, n), want in expected.items():
        vectors = list(product(range(p), repeat=n))
        triples = list(combinations(range(n), 3))
        distribution = Counter()
        for coefficients in product(range(p), repeat=len(triples)):
            t = tensor(p, n, coefficients)
            total = 0
            for x in vectors:
                b = [[sum(x[i] * t[i][j][k] for i in range(n)) % p
                      for k in range(n)] for j in range(n)]
                total += p ** (n - rank_mod(b, p))
            distribution[total] += 1
        assert dict(distribution) == want, (p, n, distribution, want)
        observed = Fraction(sum(r * count for r, count in distribution.items()), sum(distribution.values()))
        first_moment = (p + 1) * p**n - p + Fraction((p**n - 1) * (p**n - p), p ** (n - 2))
        assert observed == first_moment
        print(f"FIDELITY p={p} n={n} distribution={dict(sorted(distribution.items()))} mean={observed}")


def trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def remainder(a, b, p):
    a = trim(list(a))
    while len(a) >= len(b):
        shift = len(a) - len(b)
        c = a[-1] * pow(b[-1], -1, p) % p
        for i, v in enumerate(b):
            a[i + shift] = (a[i + shift] - c * v) % p
        trim(a)
    return a


def mulpoly(a, b, modulus, p):
    c = [0] * max(0, len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] = (c[i + j] + x * y) % p
    return remainder(c, modulus, p)


def powpoly(a, e, modulus, p):
    r = [1]
    while e:
        if e % 2:
            r = mulpoly(r, a, modulus, p)
        a = mulpoly(a, a, modulus, p)
        e //= 2
    return r


def gcdpoly(a, b, p):
    while b:
        a, b = b, remainder(a, b, p)
    return a


def irreducible(f, p):
    n = len(f) - 1
    x = remainder([0, 1], f, p)
    # Rabin criterion, using every proper divisor n/l for a prime l | n.
    prime_divisors = [l for l in range(2, n + 1)
                      if n % l == 0 and all(l % d for d in range(2, l))]
    for l in prime_divisors:
        g = powpoly(x, p ** (n // l), f, p)
        g += [0] * (max(len(g), len(x)) - len(g))
        for i, a in enumerate(x):
            g[i] = (g[i] - a) % p
        if len(gcdpoly(f, trim(g), p)) != 1:
            return False
    return powpoly(x, p**n, f, p) == x


class Field:
    def __init__(self, p, n):
        self.p, self.n, self.q = p, n, p**n
        if n == 1:
            self.modulus = [0, 1]
        else:
            self.modulus = next(list(c) + [1] for c in product(range(p), repeat=n)
                                if c[0] and irreducible(list(c) + [1], p))
        assert irreducible(self.modulus, p)
        self.coords = [tuple((a // p**i) % p for i in range(n)) for a in range(self.q)]
        self.add = [[self.encode(tuple((x + y) % p for x, y in zip(a, b)))
                     for b in self.coords] for a in self.coords]
        self.neg = [self.encode(tuple(-x % p for x in a)) for a in self.coords]
        self.mul = [[self.encode(mulpoly(a, b, self.modulus, p))
                     for b in self.coords] for a in self.coords]
        self.frob = [self.power(a, p) for a in range(self.q)]
        self.frob2 = [self.frob[self.frob[a]] for a in range(self.q)]
        self.trace = []
        for a in range(self.q):
            value, b = 0, a
            for _ in range(n):
                value = self.add[value][b]
                b = self.frob[b]
            assert b == a and value < p
            self.trace.append(value)
        self.basis = [p**i for i in range(n)]

    def encode(self, a):
        return sum(x * self.p**i for i, x in enumerate(a))

    def power(self, a, e):
        r = 1
        while e:
            if e % 2:
                r = self.mul[r][a]
            a = self.mul[a][a]
            e //= 2
        return r

    def subtract(self, a, b):
        return self.add[a][self.neg[b]]

    def abc(self, x, y):
        m, f, ff = self.mul, self.frob, self.frob2
        return (self.subtract(m[x][f[y]], m[f[x]][y]),
                self.subtract(m[ff[x]][y], m[x][ff[y]]),
                self.subtract(m[f[x]][ff[y]], m[ff[x]][f[y]]))

    def tau(self, x, y, z):
        # Direct signed determinant, with no adjoint or kernel polynomial.
        rows = [(x, y, z), (self.frob[x], self.frob[y], self.frob[z]),
                (self.frob2[x], self.frob2[y], self.frob2[z])]
        det = 0
        for perm, sg in PERMS:
            term = self.mul[self.mul[rows[0][perm[0]]][rows[1][perm[1]]]][rows[2][perm[2]]]
            det = self.add[det][term if sg == 1 else self.neg[term]]
        return self.trace[det]

    def polynomial(self, x, y):
        p, m = self.p, self.mul
        xp, xp3, xp4 = self.frob[x], self.power(x, p**3), self.power(x, p**4)
        terms = [m[xp3][self.power(y, p**4)],
                 self.neg[m[self.add[xp][xp4]][self.power(y, p**3)]],
                 m[self.add[x][xp3]][self.frob[y]], self.neg[m[xp][y]]]
        result = 0
        for term in terms:
            result = self.add[result][term]
        return result


def small_cocycle(field):
    p, n, q = field.p, field.n, field.q
    triples = list(combinations(range(n), 3))
    coeff = [field.tau(*(field.basis[i] for i in t)) for t in triples]

    def f(a, b, c):
        av, bv, cv = field.coords[a], field.coords[b], field.coords[c]
        return sum(t * av[i] * bv[j] * cv[k] for (i, j, k), t in zip(triples, coeff)) % p

    # On all basis triples, alternation of the coordinate cochain is the form.
    for x, y, z in product(field.basis, repeat=3):
        v = (x, y, z)
        alternating = sum(sg * f(*(v[i] for i in perm)) for perm, sg in PERMS) % p
        assert alternating == field.tau(x, y, z)
    checked = 0
    if q <= 9:
        for a, b, c, d in product(range(q), repeat=4):
            assert (f(b, c, d) - f(field.add[a][b], c, d)
                    + f(a, field.add[b][c], d) - f(a, b, field.add[c][d]) + f(a, b, c)) % p == 0
            checked += 1
        for a, b in product(range(q), repeat=2):
            assert f(0, a, b) == f(a, 0, b) == f(a, b, 0) == 0
    if q <= 16:
        for x, y, z in product(range(q), repeat=3):
            v = (x, y, z)
            alternating = sum(sg * f(*(v[i] for i in perm)) for perm, sg in PERMS) % p
            assert alternating == field.tau(x, y, z)
    return checked


def construction(p, n):
    f = Field(p, n)
    q = f.q
    # Nondegenerate trace pairing checked on every nonzero first argument.
    assert all(any(f.trace[f.mul[a][b]] for b in range(q)) for a in range(1, q))
    # Verify Frobenius linearity on every pair and scalar, hence all determinant
    # ingredients are checked independently of the proposed polynomial.
    for a, b in product(range(q), repeat=2):
        assert f.frob[f.add[a][b]] == f.add[f.frob[a]][f.frob[b]]
    for a in range(q):
        for c in range(p):
            assert f.frob[f.mul[c][a]] == f.mul[c][f.frob[a]]
    radical_sizes = []
    for x in range(q):
        matrix = [[f.tau(x, y, z) for z in f.basis] for y in f.basis]
        mrank = rank_mod(matrix, p)
        assert mrank % 2 == 0
        count = 0
        for y in range(q):
            readings = [f.tau(x, y, z) for z in f.basis]
            # Cross-check matrix evaluation against a direct determinant.
            assert readings == [sum(f.coords[y][j] * matrix[j][k] for j in range(n)) % p
                                for k in range(n)]
            in_radical = not any(readings)
            count += in_radical
            assert f.tau(x, x, y) == f.tau(x, y, x) == f.tau(y, x, x) == 0
            if n >= 2:
                a, b, c = f.abc(x, y)
                adjoint = f.add[f.add[f.power(a, p ** (n - 2))][f.power(b, p ** (n - 1))]][c]
                assert in_radical == (adjoint == 0)
                for z, value in zip(f.basis, readings):
                    assert value == f.trace[f.mul[adjoint][z]]
                assert f.power(adjoint, p**2) == f.polynomial(x, y)
                assert in_radical == (f.polynomial(x, y) == 0)
        assert count == p ** (n - mrank)
        assert x == 0 or count <= p**4
        radical_sizes.append(count)
    rank = sum(radical_sizes)
    assert radical_sizes[0] == q
    assert (p + 1) * q - p <= rank <= q + (q - 1) * p**4 < (p**4 + 1) * q
    if n >= 5:
        assert all(s < q for s in radical_sizes[1:])
    cocycle_checks = small_cocycle(f)
    print(f"CONSTRUCTION p={p} n={n} modulus={f.modulus} rank={rank} "
          f"nonzero_radicals={dict(sorted(Counter(radical_sizes[1:]).items()))} "
          f"pairs={q*q} cocycle_quadruples={cocycle_checks}")


def formal_exponent_identity():
    # Ordinary bivariate polynomials over F_p, with explicit sparse exponents.
    # This compares A+B^p+C^(p²) to the displayed signed P without a field
    # quotient, including characteristic 2 and characteristic 3.
    for p in (2, 3, 5, 7, 11, 13):
        a = {(1, p): 1, (p, 1): -1}
        b = {(p**2, 1): 1, (1, p**2): -1}
        c = {(p, p**2): 1, (p**2, p): -1}
        lhs = Counter()
        for polynomial, exponent in ((a, 1), (b, p), (c, p**2)):
            for (i, j), coefficient in polynomial.items():
                lhs[i * exponent, j * exponent] += pow(coefficient, exponent)
        lhs = {e: v % p for e, v in lhs.items() if v % p}
        rhs = {(p**3, p**4): 1, (p, p**3): -1, (p**4, p**3): -1,
               (1, p): 1, (p**3, p): 1, (p, 1): -1}
        rhs = {e: v % p for e, v in rhs.items() if v % p}
        assert lhs == rhs
        print(f"SIGNED_POLYNOMIAL p={p} highest_y_degree={p**4} coefficient=x^{p**3}")


def main():
    fidelity()
    formal_exponent_identity()
    cases = [(2, n) for n in range(1, 9)] + [(3, n) for n in range(1, 6)]
    cases += [(p, n) for p in (5, 7) for n in range(1, 4)]
    for p, n in cases:
        construction(p, n)
    print(f"COVERAGE fields={len(cases)} exhaustive_field_pairs={sum(p**(2*n) for p,n in cases)}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
