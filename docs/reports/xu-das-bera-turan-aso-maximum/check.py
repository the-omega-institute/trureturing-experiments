#!/usr/bin/env python3
"""Exact checks for the structural Turan-ASO conjecture (not its printed formula).

Only integers, fractions, and exact symbolic algebra are used. Radical comparisons
have rational enclosing intervals obtained from integer square roots. Enumeration
supports the preregistration; it does not prove the universally quantified claim.
"""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from math import isqrt
import sympy as S


@lru_cache(None)
def radical(q):
    """Return rational coefficient and squarefree d for sqrt(q) = coefficient*sqrt(d)."""
    q = Q(q)
    assert q > 0
    square, d = 1, 1
    for p, e in S.factorint(q.numerator * q.denominator).items():
        square *= int(p) ** (int(e) // 2)
        if e % 2:
            d *= int(p)
    return Q(square, q.denominator), d


def weighted_sum(terms):
    answer = defaultdict(Q)
    for coefficient, radicand in terms:
        c, d = radical(radicand)
        answer[d] += coefficient * c
    return {d: c for d, c in answer.items() if c}


def subtract(a, b):
    out = defaultdict(Q, a)
    for d, c in b.items():
        out[d] -= c
    return {d: c for d, c in out.items() if c}


def sign(a):
    """Certify a sign using rational intervals, without evaluating a float."""
    if not a:
        return 0
    for bits in (16, 32, 64, 128, 256, 512):
        scale = 1 << bits
        lower = upper = Q(0)
        for d, c in a.items():
            q = isqrt(d * scale * scale)
            lo = Q(q, scale)
            hi = lo if q * q == d * scale * scale else Q(q + 1, scale)
            lower += c * (lo if c > 0 else hi)
            upper += c * (hi if c > 0 else lo)
        if lower > 0:
            return 1
        if upper < 0:
            return -1
    raise AssertionError(('Sign unresolved by exact intervals', a))


def partitions(n, k, minimum=1):
    if k == 0:
        if n == 0:
            yield ()
        return
    for a in range(minimum, n // k + 1):
        for tail in partitions(n - a, k - 1, a):
            yield (a,) + tail


def balanced(n, k):
    q, r = divmod(n, k)
    return (q,) * (k - r) + (q + 1,) * r


def objective(parts):
    n = sum(parts)
    return weighted_sum((a * b, Q((n-a)**2 + (n-b)**2, 2*n-a-b-2))
                        for a, b in combinations(parts, 2))


def graph_value(n, edges):
    degrees = [0] * n
    for a, b in edges:
        degrees[a] += 1
        degrees[b] += 1
    if any(degrees[a] == degrees[b] == 1 for a, b in edges):
        return None, degrees
    return weighted_sum((1, Q(degrees[a]**2 + degrees[b]**2,
                              degrees[a] + degrees[b] - 2))
                        for a, b in edges), degrees


def colorable(n, edges, k):
    neighbors = [set() for _ in range(n)]
    for a, b in edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
    colors = [-1] * n
    def go(v, used):
        if v == n:
            return True
        for c in range(min(k, used + 1)):
            if all(colors[w] != c for w in neighbors[v]):
                colors[v] = c
                if go(v + 1, max(used, c + 1)):
                    return True
        colors[v] = -1
        return False
    return go(0, 0)


def is_turan(n, edges, k):
    """Recognize complete multipartite graphs via nonadjacency classes."""
    edges = {frozenset(e) for e in edges}
    cls = {frozenset(w for w in range(n) if w == v or frozenset((v, w)) not in edges)
           for v in range(n)}
    if len(cls) != k or sum(map(len, cls)) != n:
        return False
    if any(a & b for a, b in combinations(cls, 2)):
        return False
    return sorted(map(len, cls)) == list(balanced(n, k))


def symbolic_checks():
    m, c, u = S.symbols('m c u', positive=True)
    v = 2*m-c-u
    D = 2*m-c
    C = 2*((m-c)**2+1)
    h = S.sqrt(((m+1-u)**2+(m+1-c)**2)/v)
    rhs = (D*C/v**3-1)/h-u*(1-C/v**2)**2/(4*h**3)
    assert S.simplify(S.diff(u*h, u, 2)-rhs) == 0
    t = S.symbols('t')
    R = (4-2*t)*((1+t)**(-2)+S.Rational(1,25))
    B = 1-S.Rational(216,25)*t+S.Rational(501,25)*t**2+6*t**3+S.Rational(8,25)*t**4
    assert S.cancel((1+t)**2*((1+2*t)-4*t*(R-1))-B) == 0
    assert S.Rational(216,25)**2-4*S.Rational(501,25) == -S.Rational(3444,625)
    assert S.expand((1+6*t+t*t)**2-(1+2*t)**2*(1+t)**3 - t*(5+19*t-13*t*t-15*t**3-4*t**4)) == 0
    assert (4-2*Q(3,4))*(1/(1+Q(3,4))**2+Q(1,25)) == Q(449,490) < 1
    print('SYMBOLIC: exact curvature identity, corrected polynomial, discriminant -3444/625, lower-bound identity, R(3/4)<=449/490')
    # The uncorrected squared estimate is false, even at feasible n=6,(a,b,c)=(4,1,1).
    t = Q(1,5)
    R = (4-2*t)*(1/(1+t)**2+Q(1,25))
    bad = (1+2*t)-4*t*(R-1)**2
    assert bad == -Q(59546,78125) < 0
    print('OUTLINE CORRECTION: squared expression at m=5,s=1 is -59546/78125; the certified formula uses (1+t)^2 times the UNSQUARED bracket')


def main():
    symbolic_checks()
    fidelity = 0
    for n in range(3, 31):
        for a in range(1, n):
            b = n-a
            edges = [(i, j) for i in range(a) for j in range(a, n)]
            actual, _ = graph_value(n, edges)
            expected = weighted_sum([(a*b, Q(a*a+b*b, n-2))])
            assert actual == expected == objective((a,b))
            diff = subtract(objective(balanced(n,2)), actual)
            assert sign(diff) == (0 if abs(a-b) <= 1 else 1)
            fidelity += 1
    print(f'FIDELITY: complete bipartite formula and the source bipartite maximum/equality checked for all {fidelity} splits with 3<=n<=30')
    count = transfers = 0
    for n in range(3, 21):
        for k in range(2, n+1):
            target = objective(balanced(n,k))
            for parts in partitions(n,k):
                assert sign(subtract(target, objective(parts))) == (0 if parts == balanced(n,k) else 1), (n,k,parts)
                count += 1
                for i,j in combinations(range(k),2):
                    if parts[j] >= parts[i]+2:
                        moved = list(parts)
                        moved[j] -= 1
                        moved[i] += 1
                        assert sign(subtract(objective(tuple(moved)), objective(parts))) == 1
                        transfers += 1
    print(f'PARTITIONS: {count} positive partitions, all 3<=n<=20 and 2<=k<=n; {transfers} strict balancing transfers; no violation of maximum or equality')
    graph_count = tests = 0
    for n in range(3,6):
        possible = list(combinations(range(n),2))
        for mask in range(1 << len(possible)):
            edges = [e for i,e in enumerate(possible) if mask & (1<<i)]
            val, _ = graph_value(n, edges)
            if val is None:
                continue
            graph_count += 1
            for k in range(2,n+1):
                if not colorable(n,edges,k):
                    continue
                relation = sign(subtract(objective(balanced(n,k)), val))
                assert relation >= 0, (n,k,edges)
                assert (relation == 0) == is_turan(n,edges,k), (n,k,edges)
                tests += 1
    print(f'GRAPHS: {graph_count} labeled simple graphs without P2 components, all 3<=n<=5; {tests} colorable (G,k) cases; inequality and isomorphism equality verified')
    print('SCOPE: finite checks and exact algebraic identities; no universal proof or Lean proof is claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
