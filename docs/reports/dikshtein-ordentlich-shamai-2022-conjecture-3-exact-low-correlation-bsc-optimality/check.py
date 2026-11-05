#!/usr/bin/env python3
"""Exact finite-channel checks for DOS Conjecture 3; no universal proof.

All arithmetic and inequality certificates use integers and Fraction. Natural
logarithms of positive rationals are represented by prime-log coefficients.
Their signs are certified by rational atanh-series enclosures. No float is used.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product


@lru_cache(None)
def factors(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def logsum(terms):
    """Exact representation of sum c*ln(q), for rational c and positive q."""
    ans = {}
    for c, q in terms:
        if not c:
            continue
        assert q > 0
        for sign, n in [(1, q.numerator), (-1, q.denominator)]:
            for prime, exponent in factors(n).items():
                ans[prime] = ans.get(prime, Q(0)) + c * sign * exponent
    return tuple(sorted((p, c) for p, c in ans.items() if c))


def sub(a, b):
    return logsum([(c, Q(p)) for p, c in a] + [(-c, Q(p)) for p, c in b])


@lru_cache(None)
def log_interval(q):
    """ln(q) enclosure; 2*sum t^(2j+1)/(2j+1), with rigorous tail."""
    q = Q(q)
    exponent = 0
    while q >= 2:
        q /= 2
        exponent += 1
    while q < 1:
        q *= 2
        exponent -= 1

    def unit_interval(x):
        t = (x - 1) / (x + 1)
        n = 36
        lo = 2 * sum((t ** (2*j+1) / (2*j+1) for j in range(n)), Q(0))
        hi = lo + 2 * t ** (2*n+1) / ((2*n+1) * (1-t*t))
        return lo, hi

    lo, hi = unit_interval(q)
    l2, u2 = unit_interval(Q(2))
    if exponent >= 0:
        return lo + exponent*l2, hi + exponent*u2
    return lo + exponent*u2, hi + exponent*l2


@lru_cache(None)
def interval(expr):
    lo = hi = Q(0)
    for p, c in expr:
        a, b = log_interval(p)
        lo += c * (a if c >= 0 else b)
        hi += c * (b if c >= 0 else a)
    return lo, hi


def nonnegative(expr):
    if not expr:
        return True
    lo, hi = interval(expr)
    assert lo > 0 or hi < 0, 'insufficient exact enclosure'
    return lo > 0


def hb(q):
    return logsum([(-t, t) for t in (q, 1-q) if t])


LN2 = logsum([(Q(1), Q(2))])


def mi(joint):
    ru = [sum(row) for row in joint]
    cv = [sum(row[j] for row in joint) for j in range(len(joint[0]))]
    return logsum([(q, q / (ru[i]*cv[j]))
                   for i, row in enumerate(joint) for j, q in enumerate(row) if q])


def rate(k):
    return mi([[q/2 for q in row] for row in k])


def dsbs(p):
    return [[((1-p) if x == y else p)/2 for y in range(2)] for x in range(2)]


def joint(p, k, l):
    source = dsbs(p)
    return [[sum(source[x][y]*k[x][u]*l[y][v]
                 for x in range(2) for y in range(2))
             for v in range(len(l[0]))] for u in range(len(k[0]))]


def bsc(a):
    return ((1-a, a), (a, 1-a))


def conv(a, b):
    return a*(1-b)+(1-a)*b


def rows(n, denominator):
    return [tuple(Q(x, denominator) for x in v)
            for v in product(range(denominator+1), repeat=n) if sum(v) == denominator]


def posterior_identity(p, k, l):
    ku = [(k[0][u]+k[1][u])/2 for u in range(len(k[0]))]
    lv = [(l[0][v]+l[1][v])/2 for v in range(len(l[0]))]
    j = joint(p, k, l)
    for u, pu in enumerate(ku):
        for v, pv in enumerate(lv):
            if pu and pv:
                s = (k[0][u]-k[1][u])/(2*pu)
                t = (l[0][v]-l[1][v])/(2*pv)
                assert j[u][v] == pu*pv*(1+(1-2*p)*s*t)
            else:
                assert j[u][v] == 0


def main():
    # Fidelity: source (11), Lemma 3 product-measure identity, and the BSC
    # lower-bound value in Proposition 5, all evaluated exactly, not by fitting.
    cases = 0
    for a, b, p in product([Q(i, 8) for i in range(5)], repeat=3):
        s = dsbs(p)
        assert sum(map(sum, s)) == 1
        assert sum(s[0]) == sum(s[1]) == Q(1, 2)
        k, l = bsc(a), bsc(b)
        posterior_identity(p, k, l)
        assert mi(joint(p, k, l)) == sub(LN2, hb(conv(conv(a, p), b)))
        assert rate(k) == sub(LN2, hb(a))
        cases += 1
    print(f'FIDELITY: {cases} exact BSC value/rate/DSBS/posterior checks')

    kernels = []
    for n in (1, 2, 3):
        rr = rows(n, 4)
        kernels.extend(product(rr, repeat=2))
    assert len(kernels) == 251
    budgets = [(Q(1, 4), Q(1, 4)), (Q(1, 8), Q(1, 3))]
    comparisons = 0
    for a, b in budgets:
        cu, cv = sub(LN2, hb(a)), sub(LN2, hb(b))
        all_k = list(dict.fromkeys(kernels + [bsc(a), bsc(b)]))
        rates = {k: rate(k) for k in all_k}
        left = [k for k in all_k if nonnegative(sub(cu, rates[k]))]
        right = [k for k in all_k if nonnegative(sub(cv, rates[k]))]
        for p in [Q(49, 100), Q(99, 200)]:
            target = mi(joint(p, bsc(a), bsc(b)))
            ties = 0
            for k, l in product(left, right):
                posterior_identity(p, k, l)
                gap = sub(target, mi(joint(p, k, l)))
                assert nonnegative(gap), (a, b, p, k, l, gap)
                ties += not gap
                comparisons += 1
            print(f'GRID: a={a} b={b} p={p}; feasible={len(left)}x{len(right)}; '
                  f'comparisons={len(left)*len(right)}; exact_ties={ties}; no violation')
    # Coefficient bounds used in the proposed remainder estimate.
    for k in range(2, 129):
        assert Q(6*k*(k-1), 2*k*(2*k-1)) <= Q(3, 2)
    for k in range(1, 129):
        assert Q(6*k*k, 2*k*(2*k+1)) <= Q(3, 2)
    for M in [Q(i, 16) for i in range(1, 257)]:
        delta = 1/(1+6*M)
        assert -Q(1, 4) + 3*M*delta/(2*(1-delta)) == 0
    print(f'SETTLEMENT CHECK: {comparisons} rational finite-kernel comparisons; '
          'Taylor coefficient checks k<=128; 256 exact threshold substitutions')
    print('SCOPE: finite grid checks support the outline; they do not prove the universal conjecture')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
