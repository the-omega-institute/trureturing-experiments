#!/usr/bin/env python3
"""Exact arithmetic checks for the BDGL Question 7.1 preregistration.

No floating point, numerical root finder, or finite discretization is used as
evidence for a topological assertion. The universal algebraic certificates are
checked as polynomial identities; the remaining topology is a proof obligation.
Requires Python 3 and SymPy.
"""

from fractions import Fraction as Q
from itertools import product

import sympy as sp


def fidelity():
    """Recompute source Example 3.3, using positive iterates and discrete X."""
    points = (0, 1)
    pairs = set(product(points, repeat=2))
    diagonal = {(x, x) for x in points}
    off_diagonal = sorted(pairs - diagonal)
    neighborhoods = [
        diagonal | {p for p, bit in zip(off_diagonal, bits) if bit}
        for bits in product((0, 1), repeat=len(off_diagonal))
    ]
    f = {0: 1, 1: 1}
    # All positive powers equal f; closures in a finite discrete space do nothing.
    assert {x: f[f[x]] for x in points} == f
    rp, backward = pairs.copy(), pairs.copy()
    for alpha in neighborhoods:
        rp &= {p for p in pairs if (f[p[0]], f[p[1]]) in alpha}
        backward &= {(f[a], f[b]) for a, b in alpha}
    assert rp == pairs and backward == {(1, 1)} and backward < rp
    print("FIDELITY Example 3.3: RP = X^2; RPMinus = {(1, 1)}")
    # Identity changes the finite source example: its RPMinus becomes diagonal.
    backward_with_identity = pairs.copy()
    for alpha in neighborhoods:
        backward_with_identity &= alpha | {(f[a], f[b]) for a, b in alpha}
    assert backward_with_identity == diagonal
    print("SEMIGROUP: positive powers reproduce source; adjoining identity changes this example")


def polynomial_certificates():
    t, x, y, a, r, z = sp.symbols("t x y a r z", real=True)
    f = x + t * x * (1 - x**2)

    def identity(name, lhs, rhs):
        assert sp.Poly(sp.expand(lhs - rhs), t, x, y, a, r, z).is_zero, name

    identity("derivative positivity", sp.diff(f, x),
             (1 - 2*t) + 3*t*(1 - x)*(1 + x))
    identity("strict increase divided difference", f.subs(x, y) - f,
             (y - x) * ((1 - 2*t) + sp.Rational(3, 2)*t *
                        ((1 - x**2) + (1 - y**2)) +
                        sp.Rational(1, 2)*t*(x - y)**2))
    identity("fixed points and drift", f - x, t*x*(1 - x)*(1 + x))
    identity("oddness", f.subs(x, -x), -f)
    identity("positive endpoint barrier", 1 - f, (1 - x)*(1 - t*x*(1 + x)))
    identity("negative endpoint barrier", 1 + f, (1 + x)*(1 + t*x*(1 - x)))
    identity("forward geometric bound", (1 - t*a)*(1 - x) - (1 - f),
             t*(1 - x)*((x - a) + x**2))
    identity("inverse geometric bound",
             f.subs(x, z) - (1 + t*(1 - r**2))*z,
             t*z*(r - z)*(r + z))
    assert [f.subs(x, v) for v in (-1, 0, 1)] == [-1, 0, 1]
    assert sp.factor(f - x) == -t*x*(x - 1)*(x + 1)
    print("FAMILY: 8 exact polynomial certificates and 3 fixed-point evaluations passed")
    print("DOMAIN: 0 < t < 1/2, -1 <= x,y <= 1; strict increase has positive term 1-2t")
    print("FORWARD: 0 < a <= x < 1 gives 1-f_t(x) <= (1-t*a)*(1-x)")
    print("BACKWARD: 0 <= z <= r < 1 gives f_t(z) >= (1+t*(1-r^2))*z")


def f_exact(t, x):
    return x + t*x*(1 - x*x)


def rational_checks():
    parameters = (Q(1, 8), Q(1, 4), Q(3, 8))
    orbit_cases = 0
    for t in parameters:
        assert 0 < t < Q(1, 2)
        for x0 in (Q(-3, 4), Q(-1, 4), Q(1, 4), Q(3, 4)):
            x = x0
            for _ in range(5):
                next_x = f_exact(t, x)
                assert -1 < next_x < 1
                assert (x < next_x) if x > 0 else (next_x < x)
                a, u = abs(x0), abs(x)
                assert 1 - abs(next_x) <= (1 - t*a)*(1 - u)
                x = next_x
            orbit_cases += 1

    # Certify root brackets using only exact rational evaluations. Existence and
    # uniqueness of the inverse are supplied by continuity and strict increase.
    brackets = 0
    for t, target in product(parameters, (Q(1, 2), Q(3, 4), Q(15, 16))):
        lo, hi = Q(0), target
        for _ in range(40):
            mid = (lo + hi) / 2
            if f_exact(t, mid) < target:
                lo = mid
            else:
                hi = mid
        assert f_exact(t, lo) <= target <= f_exact(t, hi)
        assert hi - lo <= Q(1, 2**40)
        assert (1 + t*(1 - target**2))*lo <= target
        brackets += 1
    print(f"RATIONAL: {orbit_cases} forward orbits (5 steps); {brackets} inverse brackets")

    # For neighborhoods of (-1,1) of radius eps, choose target pair (-r,r).
    # Inverse contraction puts its Nth preimages in (-eps/2,eps/2).
    # These exact inequalities check the witness bounds, without approximating
    # inverse roots or asserting that finitely many epsilons prove the relation.
    bounds = []
    for t, k in product(parameters, range(1, 5)):
        eps = Q(1, 2**k)
        r = 1 - eps/2
        q = 1 / (1 + t*(1 - r*r))
        assert 0 < q < 1 and 1 - r < eps
        n, q_power = 1, q
        while r*q_power >= eps/2:
            n *= 2
            q_power *= q_power
        assert 2*r*q_power < eps
        bounds.append(n)
    print(f"WITNESS BOUNDS: {len(bounds)} exact geometric bounds; largest N = {max(bounds)}")


def reversal_and_square():
    def same_half(u, v):
        return (u <= 0 and v <= 0) or (0 <= u and 0 <= v)

    def forward(p, q):
        return same_half(p[0], q[0])

    def backward(p, q):
        return same_half(p[1], q[1])

    assert not same_half(-1, 1)
    assert backward((-1, 0), (1, 0)) and not forward((-1, 0), (1, 0))
    assert forward((0, -1), (0, 1)) and not backward((0, -1), (0, 1))
    print("RELATION TABLE: R is proper in I^2; product formulas give both incomparability witnesses")
    print("LIMIT: polynomial identities and rational bounds checked; RP identities require topology")


if __name__ == "__main__":
    fidelity()
    polynomial_certificates()
    rational_checks()
    reversal_and_square()
    print("ALL_CHECKS_OK")
