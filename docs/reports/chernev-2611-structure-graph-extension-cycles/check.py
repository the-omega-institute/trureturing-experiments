#!/usr/bin/env python3
"""Exact checks for Chernev--Chernev, arXiv:2611.07521v1, Hypothesis 6.1.

Independent standard-library implementation; no companion-repository code copied.
Field elements are base-p coefficient encodings modulo the displayed monic polynomial.
These finite checks support a prospective uniform argument, not a proof of that argument.
"""
from collections import Counter
from itertools import product


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def remainder(a, h, p):
    """Coefficient-list monic polynomial long division, independently of shift multiplication."""
    a = trim(a)
    while len(a) >= len(h):
        c, k = a[-1], len(a) - len(h)
        for j, b in enumerate(h):
            a[k + j] = (a[k + j] - c * b) % p
        a = trim(a)
    return a


def irreducible(h, p):
    # Every reducible monic polynomial has a monic divisor of degree at most half its degree.
    n = len(h) - 1
    return all(remainder(h, list(c) + [1], p)
               for d in range(1, n // 2 + 1) for c in product(range(p), repeat=d))


class Field:
    def __init__(self, p, h):
        assert p >= 2 and all(p % d for d in range(2, p))
        assert h[-1] == 1 and all(0 <= c < p for c in h)
        assert irreducible(h, p), (p, h)
        self.p, self.h, self.n = p, h, len(h) - 1
        self.q = p ** self.n
        self.coeffs = [self.decode(a) for a in range(self.q)]
        self.adds = [[self.encode([(x + y) % p for x, y in zip(a, b)])
                      for b in self.coeffs] for a in self.coeffs]
        self.negs = [self.encode([(-x) % p for x in a]) for a in self.coeffs]
        self.muls = [[self.shift_mul(a, b) for b in self.coeffs] for a in self.coeffs]
        for a, b in product(range(self.q), repeat=2):
            assert self.muls[a][b] == self.convolution_mul(a, b), (p, h, a, b)
        for a in range(1, self.q):
            assert any(self.mul(a, b) == 1 for b in range(1, self.q))

    def decode(self, a):
        c = []
        for _ in range(self.n):
            c.append(a % self.p)
            a //= self.p
        assert a == 0
        return c

    def encode(self, c):
        return sum(x * self.p ** j for j, x in enumerate(c))

    def shift_mul(self, a, b):
        # Repeated multiplication by t, reducing its top coefficient immediately.
        acc, shifted = [0] * self.n, list(a)
        for c in b:
            acc = [(x + c * y) % self.p for x, y in zip(acc, shifted)]
            lead = shifted[-1]
            shifted = [0] + shifted[:-1]
            shifted = [(x - lead * self.h[j]) % self.p for j, x in enumerate(shifted)]
        return self.encode(acc)

    def convolution_mul(self, a, b):
        c = [0] * (2 * self.n - 1)
        for i, x in enumerate(self.coeffs[a]):
            for j, y in enumerate(self.coeffs[b]):
                c[i + j] = (c[i + j] + x * y) % self.p
        return self.encode(remainder(c, self.h, self.p))

    def add(self, a, b):
        return self.adds[a][b]

    def mul(self, a, b):
        return self.muls[a][b]

    def neg(self, a):
        return self.negs[a]

    def g(self, x):
        return self.add(x, self.neg(self.mul(x, x)))

    def phi(self, v):
        a, b = v
        return self.add(a, b), self.mul(a, b)


def cycles(vertices, f):
    done, found = set(), []
    for start in vertices:
        if start in done:
            continue
        path, position, x = [], {}, start
        while x not in done and x not in position:
            position[x] = len(path)
            path.append(x)
            x = f(x)
        if x in position:
            found.append(path[position[x]:])
        done.update(path)
    return found


def verify_cycle(f, vs):
    assert len(vs) > 1 and len(set(vs)) == len(vs)
    assert all(f(v) == vs[(i + 1) % len(vs)] for i, v in enumerate(vs))


def main():
    specifications = [
        (2, [0, 1]), (3, [0, 1]), (5, [0, 1]), (7, [0, 1]),
        (2, [1, 1, 1]), (2, [1, 1, 0, 1]), (3, [1, 0, 1]),
        (2, [1, 1, 0, 0, 1]), (5, [2, 0, 1]), (3, [1, 2, 0, 1]),
        (2, [1, 0, 1, 0, 0, 1]), (7, [1, 0, 1]), (2, [1, 1, 0, 0, 0, 0, 1]),
    ]
    fields = {}
    for p, h in specifications:
        f = Field(p, h)
        fields[f.q] = f
        minus_one = f.neg(1)
        for x in range(f.q):
            assert f.phi(f.phi((minus_one, x))) == (minus_one, f.g(x))
            # Check the iterated identity through eight two-step iterates for every x.
            v, y = (minus_one, x), x
            for _ in range(8):
                v, y = f.phi(f.phi(v)), f.g(y)
                assert v == (minus_one, y)
        cs = cycles(product(range(f.q), repeat=2), f.phi)
        lengths = Counter(map(len, cs))
        assert lengths[1] == f.q
        for c in cs:
            if len(c) > 1:
                verify_cycle(f.phi, c)
        if f.n > 1 and f.q != 4:
            assert any(len(c) > 1 for c in cs)
        if f.q in (2, 3, 4, 7):
            assert len(cs) == f.q  # Source's small-ring and exceptional-field claims.
        if f.n > 1 and f.n % 2:
            s = set(range(f.p, f.q))  # Constant coefficients encode the prime field.
            assert s and all(f.g(x) in s for x in s)
            c = cycles(sorted(s), f.g)[0]
            x = c[0]
            assert x != 0
            v = (minus_one, x)
            for _ in range(2 * len(c)):
                v = f.phi(v)
            assert v == (minus_one, x) and f.phi(v) != v
        if f.p % 2 and f.n % 2 == 0:
            roots = [i for i in range(f.q) if f.mul(i, i) == minus_one]
            assert roots
            for i in roots:
                one_plus_i, one_minus_i = f.add(1, i), f.add(1, f.neg(i))
                c = [(minus_one, one_plus_i), (i, f.neg(one_plus_i)),
                     (minus_one, one_minus_i), (f.neg(i), f.add(minus_one, i))]
                verify_cycle(f.phi, c)
        # Every nonzero periodic g-point lifts; no claim that its Phi-period is exactly 2m.
        for c in cycles(range(f.q), f.g):
            for x in c:
                if x:
                    v = (minus_one, x)
                    for _ in range(2 * len(c)):
                        v = f.phi(v)
                    assert v == (minus_one, x) and f.phi(v) != v
        print(f'GF({f.q}): p={p}, n={f.n}, modulus={h}, '
              f'components={len(cs)}, cycle_lengths={dict(sorted(lengths.items()))}; '
              'two multiplications, lift and applicable invariant/four-cycle checks OK')
    # Fidelity: Example 5.2, verbatim vertex coordinates from the source.
    verify_cycle(fields[5].phi, [(4, 3), (2, 2), (4, 4), (3, 1)])
    assert fields[5].g(3) == 4 and fields[5].g(4) == 3
    print('FIDELITY: source Example 5.2 Z_5 four-cycle and g two-cycle reproduced')
    assert len(cycles(product(range(9), repeat=2), fields[9].phi)) > 9
    assert len(cycles(product(range(49), repeat=2), fields[49].phi)) > 49
    print('FIDELITY: Section 8 GF(9) and GF(49) component inequalities reproduced')
    seed = [(14, 11), (5, 8), (13, 14), (3, 10), (9, 13), (4, 15), (11, 9), (2, 12)]
    verify_cycle(fields[16].phi, seed)
    print('GF(16) eight-cycle:', seed)
    print('GF(16) exact successors:', [fields[16].phi(v) for v in seed])
    # Arithmetic case split, independent of bounded-field graph enumeration.
    for n in range(3, 257):
        d, s = n, 0
        while d % 2 == 0:
            d //= 2
            s += 1
        assert n == 2 ** s * d and d % 2 == 1
        assert (d > 1 and n % d == 0) or (d == 1 and n % 4 == 0)
    print('Characteristic-two degree split checked for 3 <= n <= 256')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
