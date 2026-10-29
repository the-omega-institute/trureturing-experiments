#!/usr/bin/env python3
"""Exact finite checks for Li--Guo, Filomat 39:13, Conjecture 2.22.

Requires Python 3 and SymPy. No floating-point arithmetic is used.
Cyclotomic polynomial remainders certify identities; rational real-root
isolation independently decides signs of the matrix spectral sums.
The finite checks do not prove an assertion about every cycle order.
"""

from fractions import Fraction
from math import ceil, floor

import sympy as sp


X, Z = sp.symbols("x z")


def cycle(n):
    assert n >= 3
    return sp.Matrix(n, n, lambda i, j: int((i - j) % n in (1, n - 1)))


class CyclotomicRing:
    """Integer coefficient arithmetic in Z[z]/(Phi_m(z))."""

    def __init__(self, m):
        self.m = m
        self.modulus = list(reversed(sp.Poly(sp.cyclotomic_poly(m, Z), Z).all_coeffs()))
        self.modulus = [int(c) for c in self.modulus]
        self.degree = len(self.modulus) - 1
        assert self.modulus[-1] == 1
        self.zero = self.scalar(0)
        self.one = self.scalar(1)
        self.z = self.reduce([0, 1])
        self.powers = [self.one]
        for _ in range(1, m):
            self.powers.append(self.mul(self.powers[-1], self.z))
        assert self.mul(self.powers[-1], self.z) == self.one

    def reduce(self, values):
        values = list(values)
        for k in range(len(values) - 1, self.degree - 1, -1):
            coefficient = values[k]
            if coefficient:
                for j, a in enumerate(self.modulus):
                    values[k - self.degree + j] -= coefficient * a
        return tuple((values + [0] * self.degree)[:self.degree])

    def scalar(self, a):
        return (int(a),) + (0,) * (self.degree - 1)

    def add(self, a, b):
        return tuple(x + y for x, y in zip(a, b))

    def neg(self, a):
        return tuple(-x for x in a)

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        values = [0] * (2 * self.degree - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    if y:
                        values[i + j] += x * y
        return self.reduce(values)

    def power(self, k):
        return self.powers[k % self.m]

    def sum(self, terms):
        value = self.zero
        for t in terms:
            value = self.add(value, t)
        return value

    def root_polynomial(self, roots):
        coefficients = [self.one]
        for root in roots:
            updated = [self.zero] * (len(coefficients) + 1)
            for k, coefficient in enumerate(coefficients):
                updated[k] = self.sub(updated[k], self.mul(root, coefficient))
                updated[k + 1] = self.add(updated[k + 1], coefficient)
            coefficients = updated
        return coefficients


def rational(value):
    return Fraction(int(sp.numer(value)), int(sp.denom(value)))


def isolated_tails(polynomial):
    """Rational enclosures of positive- and negative-eigenvalue square sums."""
    positive = [Fraction(0), Fraction(0)]
    negative = [Fraction(0), Fraction(0)]
    counts = [0, 0, 0]
    for (a, b), multiplicity in polynomial.intervals(eps=sp.Rational(1, 10**12)):
        a, b = rational(a), rational(b)
        if a == b == 0:
            counts[2] += multiplicity
        elif a >= 0 and b > 0:
            positive[0] += multiplicity * a**2
            positive[1] += multiplicity * b**2
            counts[0] += multiplicity
        elif b <= 0 and a < 0:
            negative[0] += multiplicity * b**2
            negative[1] += multiplicity * a**2
            counts[1] += multiplicity
        else:
            raise AssertionError((a, b, "root sign not isolated"))
    assert sum(counts) == polynomial.degree()
    return positive, negative, counts


def outward_gap_interval(interval, n):
    scale = 10**8
    return (Fraction(floor((interval[0] - n) * scale), scale),
            Fraction(ceil((interval[1] - n) * scale), scale))


def check_cycle(n):
    a = cycle(n)
    degrees = [sum(a.row(i)) for i in range(n)]
    assert degrees == [2] * n
    d = sp.diag(*degrees)
    q, lap = d + a, d - a
    rhs = (1 - sp.Rational(1, max(degrees))) * sum(v * (v - 1) for v in degrees)
    assert rhs == n
    assert sp.trace(a * a) == 2 * n
    assert sp.trace(q * q) == 6 * n
    characteristic = a.charpoly(X).as_poly()
    # Direct matrix polynomials, independent of the proposed trigonometric formulas.
    assert q.charpoly(X).as_expr() == sp.expand(characteristic.as_expr().subs(X, X - 2))
    assert lap.charpoly(X).as_expr() == sp.expand((-1)**n * characteristic.as_expr().subs(X, 2 - X))

    ring = CyclotomicRing(2 * n)
    # At z = exp(pi*i/n), these are exactly 2*cos(2*pi*j/n).
    roots = [ring.add(ring.power(2 * j), ring.power(-2 * j)) for j in range(n)]
    assert ring.root_polynomial(roots) == [
        ring.scalar(c) for c in reversed(characteristic.all_coeffs())
    ]
    # Thus the source's Lemma 2.5(4), 2 - 2*cos(2*pi*j/n),
    # agrees with the independently computed Laplacian characteristic polynomial.
    positive_indices = [j for j in range(n) if 4 * j < n or 4 * j > 3 * n]
    negative_indices = [j for j in range(n) if n < 4 * j < 3 * n]
    tq = ring.sum(ring.mul(roots[j], roots[j]) for j in positive_indices)
    tl = ring.sum(ring.mul(roots[j], roots[j]) for j in negative_indices)
    assert ring.add(tq, tl) == ring.scalar(2 * n)
    twice_cos = ring.add(ring.power(1), ring.power(-1))
    if n % 2 == 0:
        assert tq == tl == ring.scalar(n)
        formula = "n"
    elif n % 4 == 3:
        assert ring.mul(ring.sub(tq, ring.scalar(n - 1)), twice_cos) == ring.scalar(2)
        formula = "n-1+sec(pi/n)"
    else:
        assert ring.mul(ring.sub(tq, ring.scalar(n + 1)), twice_cos) == ring.scalar(-2)
        formula = "n+1-sec(pi/n)"

    # Independent sign decisions from the actual integer matrix polynomial.
    positive, negative, counts = isolated_tails(characteristic)
    assert counts[0] == len(positive_indices)
    assert counts[1] == len(negative_indices)
    if n % 2 == 0:
        assert positive[0] <= n <= positive[1]
        assert negative[0] <= n <= negative[1]
        verdict = "equality"
    elif n % 4 == 3:
        assert positive[0] > n
        assert negative[1] < n  # Source Theorem 2.14: strict ordinary-Laplacian inequality.
        verdict = "REFUTATION"
    else:
        assert positive[1] < n
        assert negative[0] > n  # Source Remark 2: the ordinary-Laplacian exception.
        verdict = "strictly_below"
    gap = outward_gap_interval(positive, n)
    print(f"C_{n}: Phi_degree={ring.degree}; inertia={tuple(counts)}; "
          f"TQ={formula}; R={rhs}; {verdict}; rational_gap_interval=[{gap[0]},{gap[1]}]")


def main():
    print(f"Python exact integers/Fraction; SymPy {sp.__version__}; no floating-point evidence")
    # Source Lemma 2.5(1): complete-graph Laplacian spectrum n^(n-1), 0.
    for n in range(2, 9):
        a = sp.ones(n) - sp.eye(n)
        lap = (n - 1) * sp.eye(n) - a
        assert lap.charpoly(X).as_expr() == sp.expand(X * (X - n)**(n - 1))
    print("FIDELITY: Lemma 2.5(1), K_n, 2 <= n <= 8: exact characteristic polynomials match")
    for n in range(3, 36):
        check_cycle(n)
    triangle = 2 * sp.eye(3) + cycle(3)
    assert triangle.eigenvals() == {sp.Integer(4): 1, sp.Integer(1): 2}
    assert sum(m * max(v - 2, 0)**2 for v, m in triangle.eigenvals().items()) == 4
    print("C_3: exact Q eigenvalues {4:1, 1:2}; TQ=4 > R=3; violation=1")
    print("FIDELITY: Lemma 2.5(4), Theorem 2.14 and Remark 2 checked for 3 <= n <= 35")
    print("FAMILY CHECK: all 9 orders 4*k+3 for 0 <= k <= 8 violate Conjecture 2.22")
    print("SCOPE: finite exact checks only; no Lean proof and no all-orders certification")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
