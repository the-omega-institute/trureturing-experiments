#!/usr/bin/env python3
"""Exact checks for Singh arXiv:2606.05032v1, Section 4.1.

Only integers, rational matrices, certified rational root intervals, and
certified integer nth roots enter the comparisons. No floating-point value is
used. These finite checks are error detection, not a proof of the conjecture.
"""
import argparse
from functools import lru_cache
from itertools import combinations, product
from random import Random

import sympy as S
from sympy import Matrix, Rational as Q
from sympy.core.intfunc import integer_nthroot

I2 = S.eye(2)
PAULI = [I2, Matrix([[0, 1], [1, 0]]),
         Matrix([[0, -S.I], [S.I, 0]]), S.diag(1, -1)]
J = S.diag(1, -1, -1, -1)
EXPONENTS = [Q(1, 10), Q(3, 10), Q(7, 10), Q(1), Q(3, 2),
             Q(2), Q(3), Q(5), Q(10)]
X = S.Symbol('x')
SCALE = 2**120


def rotation(q):
    """SO(3) rotation of conjugation by a quaternion divided by its norm."""
    a, b, c, d = map(S.Integer, q)
    n = a*a+b*b+c*c+d*d
    assert n > 0
    R = Matrix([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])/n
    assert R.T*R == S.eye(3) and R.det() == 1
    return R


def pauli_eigenvalues(p):
    """Convex mixture of the four Pauli conjugations, hence a unital CPTP map."""
    n = sum(p)
    assert n > 0 and all(x >= 0 for x in p)
    p0, px, py, pz = map(lambda x: Q(x, n), p)
    return [p0+px-py-pz, p0-px+py-pz, p0-px-py+pz]


def channel(lam, rin=None, rout=None):
    T = S.diag(*lam)
    if rin is not None:
        T = rout*T*rin
    return S.diag(1, T)


def reference(k, v):
    """K positive definite; sigma=K^2/Tr(K^2) has an exact inverse square root.

    sigma^(-1/2) = sqrt(h) K^(-1), h=Tr(K^2). Thus W is rational even
    though its two congruence factors need not have rational entries.
    """
    assert k > 0 and k*k > sum(a*a for a in v)
    K = k*I2 + sum((a*P for a, P in zip(v, PAULI[1:])), S.zeros(2))
    h = S.expand(S.trace(K*K))
    sigma = (K*K/h).applyfunc(S.expand)
    Ki = K.inv().applyfunc(S.simplify)
    c = h/K.det()
    assert S.trace(sigma) == 1 and sigma.det() == 1/c**2
    W = Matrix(4, 4, lambda a, b:
               S.simplify(S.trace(PAULI[a]*(h*Ki*PAULI[b]*Ki))/2))
    assert all(x.is_Rational for x in W)
    assert W == W.T
    assert J*W*J == c**2*W.inv()
    return sigma, W, c


@lru_cache(maxsize=32768)
def power_interval(x, exponent):
    """Enclose x^(m/n) between rationals with certified integer arithmetic."""
    assert x >= 0 and exponent > 0
    if x == 0:
        return Q(0), Q(0)
    m, n = int(exponent.p), int(exponent.q)
    y = x**m
    target = int(y.p)*SCALE**n
    den = int(y.q)
    k, exact = integer_nthroot(target//den, n)
    assert k**n*den <= target < (k+1)**n*den
    if k**n*den == target:
        return Q(k, SCALE), Q(k, SCALE)
    return Q(k, SCALE), Q(k+1, SCALE)


def eigen_intervals(C):
    poly = C.charpoly(X).as_poly()
    intervals = poly.intervals(eps=Q(1, 10**44))
    result = []
    for (lo, hi), multiplicity in intervals:
        assert hi >= 0
        if hi == 0:
            lo = hi = Q(0)
        assert lo >= 0
        result.extend([(lo, hi)]*multiplicity)
    assert len(result) == 4
    return list(reversed(result)), poly


def check_case(N, lam, W, c, counts):
    C = N.T*W*N
    G = c*N.T*N
    assert C == C.T
    assert N*J == J*N
    # Fidelity: the source's second-moment comparison, (4.7)-(4.10).
    assert S.trace(C) >= 2*S.trace(N.T*N)
    # Fidelity: I/2 gives s_a(N,I/2)^2 = 2 s_a(N)^2.
    assert N.T*(2*S.eye(4))*N == 2*(N.T*N)
    if N.det() != 0:
        assert G*C.inv()*G == J*C*J
        counts['riccati'] += 1
    counts['cases'] += 1
    if W == 2*S.eye(4):
        assert c == 2 and C == G
        counts['equal_cases'] += 1
        counts['powers'] += len(EXPONENTS)
        counts['prefixes'] += 4
        return
    roots, poly = eigen_intervals(C)
    target = sorted([c] + [c*a*a for a in lam], reverse=True)
    rank = sum(x != 0 for x in target)
    assert sum(hi != 0 for lo, hi in roots) == rank
    for k in range(1, 5):
        rhs = S.prod(target[:k])
        if k > rank:
            assert rhs == 0
        elif k == rank:
            # Product of all nonzero roots is exactly a characteristic coefficient.
            pseudodet = (-1)**rank * poly.nth(4-rank)
            assert pseudodet >= rhs
        else:
            if S.prod(lo for lo, hi in roots[:k]) < rhs:
                # An exact equality can straddle any nonzero-width root interval.
                # A rational eigenvalue of the kth compound is at most its
                # largest eigenvalue, the product of the k largest roots of C.
                assert S.prod(hi for lo, hi in roots[:k]) >= rhs
                subsets = list(combinations(range(4), k))
                compound = Matrix(len(subsets), len(subsets),
                                  lambda a, b: C.extract(subsets[a], subsets[b]).det())
                assert (compound-rhs*S.eye(len(subsets))).det() == 0
                counts['exact_compound_fallbacks'] += 1
        counts['prefixes'] += 1
    for t in EXPONENTS:
        q = t/2
        lower = sum(power_interval(lo, q)[0] for lo, hi in roots)
        upper = sum(power_interval(x, q)[1] for x in target)
        # Tr(C^(t/2)) >= Tr((c N*N)^(t/2)), a stronger comparison
        # than the conjecture's comparison to the maximally mixed state.
        assert lower >= upper, (counts['cases'], t, lower-upper)
        counts['powers'] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--random-cases', type=int, default=1200)
    args = parser.parse_args()
    rng = Random(439)
    counts = dict(cases=0, powers=0, prefixes=0, riccati=0, equal_cases=0,
                  exact_compound_fallbacks=0)
    # Complete small Pauli probability grids: every nonnegative integer
    # 4-tuple with total 1,...,4, at three Bloch axes and three polarizations.
    grid = 0
    for n in range(1, 5):
        for p in product(range(n+1), repeat=4):
            if sum(p) != n:
                continue
            lam = pauli_eigenvalues(p)
            N = channel(lam)
            for axis in range(3):
                for k, magnitude in [(2, 1), (10, 9), (100, 99)]:
                    v = [0, 0, 0]
                    v[axis] = magnitude
                    sigma, W, c = reference(k, v)
                    check_case(N, lam, W, c, counts)
                    grid += 1
    # Rational random unitary pre/post conjugations of Pauli channels.
    # Forced faces and deficient ranks occur in addition to generic channels.
    rank_counts = {}
    for index in range(args.random_cases):
        p = [rng.randrange(1, 30) for _ in range(4)]
        if index % 10 == 0:
            p = [1, 1, 1, 1]  # rank 1
        elif index % 10 == 1:
            p = [1, 1, 0, 0]  # rank 2
        elif index % 10 == 2:
            p = [2, 1, 1, 0]  # rank 3
        elif index % 10 == 3:
            p = [1, 0, 0, 0]  # unitary extreme point
        lam = pauli_eigenvalues(p)
        rank = 1+sum(a != 0 for a in lam)
        rank_counts[rank] = rank_counts.get(rank, 0)+1
        qin = [1]+[rng.randrange(-5, 6) for _ in range(3)]
        qout = [1]+[rng.randrange(-5, 6) for _ in range(3)]
        N = channel(lam, rotation(qin), rotation(qout))
        v = [rng.randrange(-5, 6) for _ in range(3)]
        k = sum(abs(a) for a in v)+1
        if index % 24 == 0:
            v = [0, 0, 0]
        sigma, W, c = reference(k, v)
        check_case(N, lam, W, c, counts)
    assert power_interval(Q(0), Q(1, 20)) == (0, 0)
    print('arXiv:2606.05032v1, Section 4.1, equations (4.7)-(4.12)')
    print('Arithmetic: rational matrices, rational Sturm root intervals, certified integer nth roots')
    print('Fidelity: mixed-state Gram normalization and source second-moment inequality passed')
    print('All-small-cases axis grid:', grid)
    print('Random CPTP cases: seed=439, count='+str(args.random_cases)+', ranks='+str(sorted(rank_counts.items())))
    print('Exponents t:', ', '.join(map(str, EXPONENTS)))
    print('Counts:', ', '.join(f'{k}={v}' for k, v in counts.items()))
    print('Every tested sigma also passed J W J = c^2 W^(-1) exactly')
    print('Finite checks only; no universal proof or Lean kernel verification is claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
