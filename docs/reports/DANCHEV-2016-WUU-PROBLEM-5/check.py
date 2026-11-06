#!/usr/bin/env python3
"""Exact finite diagnostics for Danchev, Weakly UU rings, Problem 5.

Only Python integers are used. Finite checks do not prove the universal claim.
"""

from functools import lru_cache
from itertools import product
from math import gcd


class Matrices:
    def __init__(self, d, modulus, triangular=False):
        self.d = d
        self.m = modulus
        self.triangular = triangular
        self.zero = (0,) * (d * d)
        self.one = tuple(int(i == j) for i in range(d) for j in range(d))
        positions = [(i, j) for i in range(d) for j in range(d)
                     if not triangular or i <= j]
        self.elements = []
        for values in product(range(modulus), repeat=len(positions)):
            a = [0] * (d * d)
            for (i, j), v in zip(positions, values):
                a[i * d + j] = v
            self.elements.append(tuple(a))

    def add(self, a, b):
        return tuple((x + y) % self.m for x, y in zip(a, b))

    def scale(self, t, a):
        return tuple(t * x % self.m for x in a)

    def sub(self, a, b):
        return self.add(a, self.scale(-1, b))

    def mul(self, a, b):
        d, m = self.d, self.m
        return tuple(sum(a[i * d + k] * b[k * d + j] for k in range(d)) % m
                     for i in range(d) for j in range(d))

    def power(self, a, n):
        out = self.one
        for _ in range(n):
            out = self.mul(out, a)
        return out

    def is_unit(self, a):
        if self.triangular:
            return all(gcd(a[i * self.d + i], self.m) == 1
                       for i in range(self.d))
        assert self.d == 2
        return gcd(a[0] * a[3] - a[1] * a[2], self.m) == 1

    def inverse(self, a):
        assert self.is_unit(a)
        if not self.triangular:
            z = pow((a[0] * a[3] - a[1] * a[2]) % self.m, -1, self.m)
            return self.scale(z, (a[3], -a[1], -a[2], a[0]))
        d = self.d
        diagonal_inverse = tuple(pow(a[i * d + i], -1, self.m) if i == j else 0
                                 for i in range(d) for j in range(d))
        n = self.sub(self.mul(diagonal_inverse, a), self.one)
        out, term = self.zero, self.one
        for _ in range(d):
            out = self.add(out, term)
            term = self.mul(term, self.scale(-1, n))
        assert term == self.zero
        return self.mul(out, diagonal_inverse)

    @lru_cache(maxsize=None)
    def nil_index(self, a):
        # Decide nilpotence by the exact finite orbit of successive powers.
        # A repeated nonzero power can never subsequently reach zero.
        out, n, seen = a, 1, set()
        while out not in seen:
            if out == self.zero:
                return n
            seen.add(out)
            out = self.mul(out, a)
            n += 1
        return None

    def idempotents(self):
        return [a for a in self.elements if self.mul(a, a) == a]

    def wuu(self):
        return all(self.nil_index(self.sub(a, self.one)) is not None or
                   self.nil_index(self.add(a, self.one)) is not None
                   for a in self.elements if self.is_unit(a))


def fidelity_and_claim_instances():
    """Example 2.2(c), plus the quotient-centrality test in its WUU cases."""
    for d, p, expected in [(2, 2, True), (3, 2, True),
                            (2, 3, False), (3, 3, False)]:
        ring = Matrices(d, p, triangular=True)
        units = [a for a in ring.elements if ring.is_unit(a)]
        for a in units:
            b = ring.inverse(a)
            assert ring.mul(a, b) == ring.one == ring.mul(b, a)
        actual = ring.wuu()
        assert actual == expected
        assert len(units) == (p - 1) ** d * p ** (d * (d - 1) // 2)
        print(f"FIDELITY T_{d}(F_{p}): elements={len(ring.elements)} "
              f"units={len(units)} WUU={actual}")
        if not actual:
            continue
        # The finite Jacobson criterion, checked against every r, independently
        # recovers the strictly upper triangular radical in these examples.
        radical = {a for a in ring.elements
                   if all(ring.is_unit(ring.sub(ring.one, ring.mul(r, a)))
                          for r in ring.elements)}
        strict_upper = {a for a in ring.elements
                        if all(a[i * d + i] == 0 for i in range(d))}
        assert radical == strict_upper
        tests = 0
        for e in ring.idempotents():
            f = ring.sub(ring.one, e)
            for r in ring.elements:
                assert ring.sub(ring.mul(e, r), ring.mul(r, e)) in radical
                assert ring.mul(ring.mul(e, r), f) in radical
                assert ring.mul(ring.mul(f, r), e) in radical
                tests += 1
        print(f"CLAIM_INSTANCE T_{d}(F_{p}): radical={len(radical)} "
              f"idempotents={len(ring.idempotents())} commutator_tests={tests}")


def opposite_corner_checks(ring, label):
    pairs = implications = nonzero_nilpotents = power_tests = 0
    for e in ring.idempotents():
        f = ring.sub(ring.one, e)
        left = {ring.mul(ring.mul(e, r), f) for r in ring.elements}
        right = {ring.mul(ring.mul(f, r), e) for r in ring.elements}
        for a in left:
            assert ring.mul(a, a) == ring.zero
            for b in right:
                assert ring.mul(b, b) == ring.zero
                ab = ring.mul(a, b)
                c = ring.add(ab, ring.mul(b, a))
                premise = True
                for t in [1, 2]:
                    ta = ring.scale(t, a)
                    u = ring.mul(ring.add(ring.one, ta), ring.add(ring.one, b))
                    v = ring.mul(ring.sub(ring.one, b), ring.sub(ring.one, ta))
                    assert ring.mul(u, v) == ring.one == ring.mul(v, u)
                    tc = ring.scale(t, c)
                    assert ring.sub(ring.add(u, v), ring.scale(2, ring.one)) == tc
                    premise = premise and (ring.nil_index(tc) is not None or
                        ring.nil_index(ring.add(tc, ring.scale(4, ring.one))) is not None)
                if premise:
                    assert ring.nil_index(c) is not None
                    assert ring.nil_index(ab) is not None
                    implications += 1
                    nonzero_nilpotents += int(ab != ring.zero)
                pairs += 1
            # This checks the separate passage from are to ar, without WUU.
            for r in ring.elements:
                ar = ring.mul(a, r)
                are = ring.mul(ar, e)
                for n in range(4):
                    assert ring.power(ar, n + 1) == ring.mul(ring.power(are, n), ar)
                    power_tests += 1
                n = ring.nil_index(are)
                if n is not None:
                    assert ring.power(ar, n + 1) == ring.zero
    print(f"CORNERS {label}: idempotents={len(ring.idempotents())} "
          f"pairs={pairs} premise_true={implications} "
          f"nonzero_nilpotent_ab={nonzero_nilpotents} power_tests={power_tests}")
    return nonzero_nilpotents


def main():
    fidelity_and_claim_instances()
    for m in range(2, 8):
        opposite_corner_checks(Matrices(2, m), f"M_2(Z/{m}Z)")
    nonzero = opposite_corner_checks(Matrices(3, 2, triangular=True), "T_3(F_2)")
    assert nonzero > 0, "Include nonzero nilpotent corner products."
    print("BOUNDARY finite instances and identities only; universal theorem not certified")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
