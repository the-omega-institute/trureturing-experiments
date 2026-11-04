#!/usr/bin/env python3
"""Exact checks for Caracciolo--Vanoni, arXiv:2610.26536v1, (6.30).

Uses only integer arithmetic and fractions.Fraction. Square-root sums are
compared by certified rational intervals, never floating-point approximations.
These finite checks support a proposed general argument, not a proof of (6.30).
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import comb, factorial, isqrt


def moment_product(n, k, p):
    """Gamma(N+1)Gamma(k+p)/(Gamma(N+p+1)Gamma(k)), exact for rational p."""
    ans = Q(factorial(n), factorial(k - 1))
    for j in range(k, n + 1):
        ans /= p + j
    return ans


def moment_integral(n, k, p):
    """Integrate t^p against the Beta(k,N+1-k) polynomial density."""
    c = Q(factorial(n), factorial(k - 1) * factorial(n - k))
    return c * sum(((-1) ** j * Q(comb(n - k, j), 1) / (p + k + j)
                    for j in range(n - k + 1)), Q(0))


def reference_edges(n):
    """Literal one-based Section 6.2 edge sum for h-bar."""
    m = n // 2
    return [(m + 1, 1), (n, m)] + [
        e for k in range(2, m + 1)
        for e in [(m + k - 1, k), (m + k, k - 1)]]


def check_reference(n):
    edges = reference_edges(n)
    assert len(edges) == n
    adj = {i: [] for i in range(1, n + 1)}
    for i, j in edges:
        assert 1 <= i <= n and 1 <= j <= n and i != j
        adj[i].append(j)
        adj[j].append(i)
    assert all(len(v) == 2 for v in adj.values())
    seen, todo = set(), [1]
    while todo:
        i = todo.pop()
        if i not in seen:
            seen.add(i)
            todo.extend(adj[i])
    assert len(seen) == n
    m = n // 2
    assert Counter(abs(i - j) for i, j in edges) == Counter({m - 1: m - 1, m: 2, m + 1: m - 1})


def radical(d):
    """Return integers a,b with sqrt(d)=a*sqrt(b), b square-free."""
    a, b, q = 1, d, 2
    while q * q <= b:
        while b % (q * q) == 0:
            a *= q
            b //= q * q
        q += 1
    return a, b


def cost(points, edges):
    c = Counter()
    for i, j in edges:
        a, b = radical(abs(points[i] - points[j]))
        c[b] += a
    return dict(c)


def lincomb(*terms):
    c = Counter()
    for a, expr in terms:
        for b, v in expr.items():
            c[b] += a * v
    return {b: v for b, v in c.items() if v}


def sign(expr):
    """Exact sign by dyadic isolation; zero is recognized symbolically."""
    if not expr:
        return 0
    for bits in (8, 16, 32, 64, 128, 256):
        scale, lo, hi = 1 << bits, 0, 0
        for b, c in expr.items():
            v = isqrt(b * scale * scale)
            w = v if v * v == b * scale * scale else v + 1
            lo += c * (v if c > 0 else w)
            hi += c * (w if c > 0 else v)
        if lo > 0:
            return 1
        if hi < 0:
            return -1
    raise ArithmeticError(('unresolved radical sign', expr))


def matchings(vertices):
    if not vertices:
        yield []
        return
    a, *rest = vertices
    for j, b in enumerate(rest):
        for tail in matchings(rest[:j] + rest[j + 1:]):
            yield [(a, b)] + tail


def check_pointwise(n, grid):
    tours = [list(zip(t, t[1:] + t[:1]))
             for tail in permutations(range(1, n))
             if tail[0] < tail[-1] for t in [(0,) + tail]]
    pairs = list(matchings(list(range(n))))
    tour_count = matching_count = samples = 0
    for points in combinations(range(grid), n):
        samples += 1
        m = n // 2
        crossing = cost(points, [(i, m + i) for i in range(m)])
        for edges in pairs:
            assert sign(lincomb((1, crossing), (-1, cost(points, edges)))) >= 0
            matching_count += 1
        ref = cost(points, [(i - 1, j - 1) for i, j in reference_edges(n)])
        optimum = ref
        for edges in tours:
            value = cost(points, edges)
            assert sign(lincomb((2, crossing), (-1, value))) >= 0
            if sign(lincomb((1, value), (-1, optimum))) > 0:
                optimum = value
            tour_count += 1
        assert sign(lincomb((1, optimum), (-1, ref))) >= 0
        assert sign(lincomb((2, crossing), (-1, optimum))) >= 0
    print(f'POINTWISE p=1/2 N={n} integer_grid=0..{grid-1} '
          f'configurations={samples} tours={tour_count} matchings={matching_count} OK')
    return tour_count, matching_count


def main():
    ps = [Q(1, 5), Q(1, 3), Q(1, 2), Q(2, 3), Q(4, 5)]
    fidelity = 0
    for n in range(3, 15):
        for k in range(1, n):
            for p in ps:
                a = moment_product(n, k, p)
                assert a == moment_integral(n, k, p)
                assert 0 < a < 1
                fidelity += 1
    print(f'FIDELITY Eq.(6.4) beta integral = gamma-product: {fidelity} exact cases OK')
    envelopes = 0
    for n in range(4, 42, 2):
        check_reference(n)
        m = n // 2
        for p in ps:
            lower, mid, upper = (moment_product(n, k, p) for k in [m - 1, m, m + 1])
            ref = (m - 1) * (lower + upper) + 2 * mid
            from_edges = sum((moment_product(n, abs(i - j), p)
                              for i, j in reference_edges(n)), Q(0))
            assert ref == from_edges
            envelope = n * mid - ref
            identity = mid * p * (1 - p) * (m - 1) / (m * (m - 1 + p))
            assert envelope == identity
            assert 0 <= envelope <= 2 * p * (1 - p) / n <= Q(1, 2 * n)
            envelopes += 1
    print('FIDELITY Eqs.(6.23)-(6.24) reference edge multiplicities and expectation: '
          f'{envelopes} exact cases OK')
    print('ENVELOPE N*mu_m-E(hbar) = '
          'mu_m*p*(1-p)*(m-1)/(m*(m-1+p)): '
          f'{envelopes} exact cases OK')
    print('BOUND envelope <= 2*p*(1-p)/N <= 1/(2*N): all envelope cases OK')
    # A symbolic rational-function identity, verified by integer polynomial coefficients:
    # m(m-1+p)[2-(m-1)/(m-1+p)-(m+p)/m] = p(1-p).
    # Expand independently using convolution on (m,p) exponents.
    def mul(a, b):
        out = Counter()
        for (i, j), x in a.items():
            for (k, l), y in b.items():
                out[i + k, j + l] += x * y
        return dict(out)
    M, P, A, B = {(1, 0): 1}, {(0, 1): 1}, {(1, 0): 1, (0, 0): -1}, {(1, 0): 1, (0, 0): -1, (0, 1): 1}
    left = Counter({k: 2*v for k, v in mul(M, B).items()})
    for term in [mul(M, A), mul({(1, 0): 1, (0, 1): 1}, B)]:
        for k, v in term.items():
            left[k] -= v
    assert {k: v for k, v in left.items() if v} == {(0, 1): 1, (0, 2): -1}
    print('SYMBOLIC neighboring-moment coefficient identity: integer polynomial expansion OK')
    totals = [check_pointwise(4, 7), check_pointwise(6, 8), check_pointwise(8, 9)]
    print(f'TOTAL exact tour comparisons={sum(t for t, _ in totals)} '
          f'matching comparisons={sum(m for _, m in totals)}')
    print('SCOPE finite checks and a symbolic identity; no integration of the optimal tour, '
          'no Lean proof, and no universal settlement claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
