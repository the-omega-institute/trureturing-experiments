#!/usr/bin/env python3
"""Exact finite checks for Dimitrov--Shapiro, arXiv:2606.17864v1.

Fidelity: Lemma 3.1's reciprocal-tail bound and Proposition 3.6's
nonmultiplication example. Endpoint: every root multiset on the declared
grid through degree six, with exact reconstruction of the proposed family.
These checks do not certify the analytic assertion about entire functions.
"""

from fractions import Fraction as F
from itertools import combinations_with_replacement


def separated(roots, q):
    """Definition 1.2, with distinct indices and unrestricted zeros at zero."""
    assert F(0) < q < F(1)
    for i, x in enumerate(roots):
        for j, y in enumerate(roots):
            if i != j and x * y > 0 and abs(x) <= abs(y):
                if abs(x) / abs(y) > q:
                    return False
    return True


def critical_ratio(roots):
    return max(
        (abs(x) / abs(y)
         for i, x in enumerate(roots)
         for j, y in enumerate(roots)
         if i != j and x * y > 0 and abs(x) <= abs(y)),
        default=F(0),
    )


def multiply(p, r):
    out = [F(0)] * (len(p) + len(r) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(r):
            out[i + j] += a * b
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(out)


def polynomial(roots, scalar):
    p = (scalar,)
    for x in roots:
        p = multiply(p, (-x, F(1)))
    return p


def endpoint_parameters(roots, scalar):
    negative = [x for x in roots if x < 0]
    positive = [x for x in roots if x > 0]
    assert len(negative) <= 1 and len(positive) <= 1
    m = roots.count(F(0))
    c = scalar
    for x in negative + positive:
        c *= -x
    a = -F(1) / negative[0] if negative else F(0)
    b = F(1) / positive[0] if positive else F(0)
    assert c != 0 and a >= 0 and b >= 0
    p = (F(0),) * m + (c,)
    return multiply(multiply(p, (F(1), a)), (F(1), -b))


def main():
    # Proposition 3.6: 1-z and 1-z/y are individually admissible, their
    # product is not; use q=1/2 and y=3/2, between 1 and q^{-1}.
    q, y = F(1, 2), F(3, 2)
    assert F(1) < y < F(1) / q
    assert separated((F(1),), q) and separated((y,), q)
    assert not separated((F(1), y), q)
    assert critical_ratio((F(1), y)) == F(2, 3)
    print("FIDELITY Proposition 3.6: q=1/2, roots 1 and 3/2; ratio=2/3>1/2")

    grid = tuple(F(i) for i in range(-4, 5))
    qs = sorted({F(n, d) for d in range(2, 9) for n in range(1, d)})
    counts = dict(multisets=0, endpoint=0, excluded=0, tail_checks=0,
                  reconstructions=0, type_i=0)
    for degree in range(7):
        for roots in combinations_with_replacement(grid, degree):
            counts["multisets"] += 1
            nminus = sum(x < 0 for x in roots)
            nplus = sum(x > 0 for x in roots)
            target = nminus <= 1 and nplus <= 1
            ratio = critical_ratio(roots)
            assert (ratio == 0) == target
            if not target:
                witness_q = ratio / 2
                assert F(0) < witness_q < F(1)
                assert not separated(roots, witness_q)
                counts["excluded"] += 1
            else:
                counts["endpoint"] += 1
                for q in qs:
                    assert separated(roots, q)
                for scalar in (F(-2, 3), F(1), F(3, 2)):
                    assert polynomial(roots, scalar) == endpoint_parameters(roots, scalar)
                    counts["reconstructions"] += 1
            one_sign = nminus == 0 or nplus == 0
            assert (target and one_sign) == (nminus + nplus <= 1)
            if target and one_sign:
                counts["type_i"] += 1
            for q in qs:
                assert separated(roots, q) == (ratio <= q)
                if separated(roots, q):
                    # Lemma 3.1: exact finite root sums, both half-axes,
                    # several radii and powers; zero roots are excluded.
                    for radius in (F(1, 2), F(1), F(2), F(3), F(5)):
                        for power in (1, 2, 3):
                            tail = sum((F(1) / abs(x) ** power
                                        for x in roots if abs(x) >= radius), F(0))
                            bound = F(2) / ((F(1) - q ** power) * radius ** power)
                            assert tail <= bound
                            counts["tail_checks"] += 1

    # Boundary and sharpness checks: opposite signs, repeated nonzero
    # entries, arbitrary sampled origin multiplicities, and the zero scalar.
    for m in range(21):
        roots = (F(0),) * m + (F(-1), F(1, 2))
        for q in (F(1, 10**6), F(1, 2), F(999999, 10**6)):
            assert separated(roots, q)
        assert polynomial(roots, F(0)) == (F(0),)
    assert not separated((F(1), F(1)), F(999999, 10**6))
    assert len(qs) == 21 and counts["multisets"] == 5005
    print("GRID: roots in {-4,-3,-2,-1,0,1,2,3,4}; degrees 0..6; 21 rational q values")
    print("COUNTS: " + ", ".join(f"{key}={value}" for key, value in counts.items()))
    print("BOUNDARIES: opposite-sign roots survive; nonzero multiplicities fail; zero scalar included")
    print("SCOPE: finite exact checks support the outline; no analytic or Lean proof is claimed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
