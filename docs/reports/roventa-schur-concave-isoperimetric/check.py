#!/usr/bin/env python3
"""Exact checks for Roventa, JIA 2012:159, Section 2's final question.

All evidence uses integers and Fraction. Root enclosures certify rational
p-th-power inequalities; no numerical root or floating-point comparison is used.
This checks finite instances and algebraic identities, not the universal theorem.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product


def mean(x):
    return sum(x, Q(0)) / len(x)


def replace(x, i):
    return tuple(mean(x) if j == i else a for j, a in enumerate(x))


def negative_moment(x, p):
    assert p > 0 and all(a > 0 for a in x)
    return mean(tuple(a ** (-p) for a in x))


@lru_cache(maxsize=None)
def root_bounds(x, p):
    """Enclose M_{-p}(x), defined by M^p = 1 / negative_moment."""
    target = 1 / negative_moment(x, p)
    lo, hi = min(x), max(x)
    assert 0 < lo <= hi and lo**p <= target <= hi**p
    for _ in range(48):
        mid = (lo + hi) / 2
        if mid**p <= target:
            lo = mid
        else:
            hi = mid
    assert lo**p <= target <= hi**p
    return lo, hi


def failure_certificate(x, p):
    n, mu = len(x), mean(x)
    flo, _ = root_bounds(x, p)
    hhi = sum((root_bounds(replace(x, i), p)[1] for i in range(n)), Q(0))
    # F(mu,...,mu)=mu, H(mu,...,mu)=n*mu exactly.
    # These rational bounds imply A^(n-1) > B^n, hence strict failure of (2.3).
    gap = (flo / mu) ** (n - 1) - (hhi / (n * mu)) ** n
    assert gap > 0, (x, p)
    return gap


def fidelity():
    count = 0
    # Theorem 5: f(t)=t/(1+t), positive and concave on the positive axis.
    # Its derivatives are 1/(1+t)^2 and -2/(1+t)^3.
    f = lambda t: t / (1 + t)
    for n in range(2, 7):
        for data in product((1, 2, 4), repeat=n):
            x = tuple(map(Q, data))
            fx = sum(map(f, x), Q(0))
            fm = f(mean(x))
            hx = sum((sum(map(f, replace(x, i)), Q(0)) for i in range(n)), Q(0))
            assert hx == (n - 1) * fx + n * fm  # Source's equality in (2.6).
            a, b = fx / (n * fm), hx / (n * n * fm)
            assert 0 < a <= 1 and a ** (n - 1) <= b**n
            count += 1
    print(f'FIDELITY: Theorem 5 and equality (2.6), {count} exact rational cases, n=2..6')


def admissibility_checks():
    count = 0
    for n in range(3, 9):
        x, p = tuple(Q(i + 1) for i in range(n)), n * n
        s = negative_moment(x, p)
        a = min(x)
        assert a**(-p) / n <= s <= a**(-p)  # Minimum envelope before reversing power.
        for theta in (Q(0), Q(1, 3), Q(1, 2), Q(2, 3), Q(1)):
            mixed = list(x)
            mixed[0] = theta * x[0] + (1 - theta) * x[-1]
            mixed[-1] = (1 - theta) * x[0] + theta * x[-1]
            assert sum(mixed) == sum(x)
            assert negative_moment(tuple(mixed), p) <= s
            count += 1
        # Exact Schur-Ostrowski sign after factoring out its positive common factor.
        for i in range(n):
            for j in range(n):
                sign = (x[i] - x[j]) * (x[i]**(-p-1) - x[j]**(-p-1))
                assert sign == 0 if i == j else sign < 0
        # Weighted Cauchy-Schwarz identity controlling the Hessian in the outline.
        v = tuple(Q((-1)**i * (i + 2)) for i in range(n))
        w = tuple(t**(-p) for t in x)
        lhs = sum(w) * sum((w[i] * (v[i]/x[i])**2 for i in range(n)), Q(0))
        lhs -= sum((w[i] * v[i]/x[i] for i in range(n)), Q(0))**2
        rhs = sum((w[i]*w[j]*(v[i]/x[i]-v[j]/x[j])**2
                   for i in range(n) for j in range(i+1, n)), Q(0))
        assert lhs == rhs and rhs >= 0
    print(f'ADMISSIBILITY: {count} rational T-transforms; strict derivative signs and Hessian identity, n=3..8')


def counterfamilies():
    for n in range(3, 17):
        x = (Q(1),) * (n-1) + (Q(n+1),)
        assert mean(x) == 2 and all(min(replace(x, i)) == 1 for i in range(n))
        for p in (n*n, 2*n*n):
            # This integer condition is equivalent to p > n*log_2(n).
            assert n**n < 2**p
            gap = failure_certificate(x, p)
            # The loose envelope alone already certifies the separation.
            assert Q(n**n, 2**(p*n)) < Q(1, 2**(p*(n-1)))
        print(f'FAMILY: n={n}, p={n*n} and {2*n*n}; exact root enclosures give A^(n-1)>B^n')
    print('FAMILY: 28 strict counterexamples, n=3..16; universal n and real p require the paper proof')
    for n in range(3, 9):
        x = tuple(Q(1) + Q(i, 2**18) for i in range(n-1)) + (Q(n+1),)
        assert len(set(x)) == n
        failure_certificate(x, n*n)
    print('ROBUSTNESS: 6 all-distinct rational inputs, n=3..8, p=n^2')
    for x in [(Q(1), Q(1), Q(3,2)), (Q(2), Q(2), Q(3), Q(5)),
              (Q(3,7), Q(3,7), Q(9,7), Q(15,7), Q(12,7))]:
        n, p = len(x), len(x)**2
        while (mean(x)/min(x))**p <= n**n:
            p *= 2
        assert all(min(replace(x, i)) == min(x) for i in range(n))
        failure_certificate(x, p)
    print('GENERAL FAMILY: 3 repeated-minimum inputs; exact threshold (mu/a)^p>n^n')


if __name__ == '__main__':
    fidelity()
    admissibility_checks()
    counterfamilies()
    print('ALL_CHECKS_OK')
