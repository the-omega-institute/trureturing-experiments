#!/usr/bin/env python3
"""Exact supplementary checks; finite searches do not prove the infinite claim.

Use only integer arithmetic and Fraction pairs (p,q) representing p+q*sqrt(D).
For t>=0, floor(t*alpha)=((1-a)*t+isqrt((a*a+4)*t*t))//4 exactly:
the discarded remainder is in [0,1) and cannot cross a multiple of four.
"""
from fractions import Fraction as Q
from math import isqrt


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def scale(c, x):
    return c*x[0], c*x[1]


def mul(x, y, D):
    return x[0]*y[0]+D*x[1]*y[1], x[0]*y[1]+x[1]*y[0]


def sign(x, D):
    p, q = x
    if not q:
        return (p > 0)-(p < 0)
    if not p:
        return (q > 0)-(q < 0)
    if (p > 0) == (q > 0):
        return 1 if p > 0 else -1
    d = p*p-D*q*q
    return ((d > 0)-(d < 0))*(1 if p > 0 else -1)


def floor_alpha(a, t):
    assert a >= 4 and t >= 0
    return ((1-a)*t+isqrt((a*a+4)*t*t))//4


def chi(a, n):
    c = floor_alpha(a, 2*n)-2*floor_alpha(a, n)
    assert c in (0, 1)
    return c


def fib_gaps(N):
    f, g, out = 1, 1, set()
    while f <= N:
        out.add(f)
        f, g = g, f+g
    return sorted(out)


def proper(N, r):
    gaps, colors = set(fib_gaps(N-1)), []
    def extend():
        if len(colors) == N:
            return tuple(colors)
        v = len(colors)+1
        forbidden = {colors[u-1] for u in range(1, v) if v-u in gaps}
        for c in range(min(r, max(colors, default=-1)+2)):
            if c not in forbidden:
                colors.append(c)
                ans = extend()
                if ans is not None:
                    return ans
                colors.pop()
        return None
    return extend()


def longest(N, gaps, coloring, r):
    """DAG recurrence checks every allowed gap, including unequal adjacent gaps."""
    colors = [0]+[coloring(n) for n in range(1, N+1)]
    assert all(0 <= c < r for c in colors[1:])
    dp, maxima, counts = [0]*(N+1), [0]*r, [0]*r
    for v in range(1, N+1):
        c, best = colors[v], 1
        for d in gaps:
            if d >= v:
                break
            if colors[v-d] == c:
                best = max(best, dp[v-d]+1)
        dp[v] = best
        maxima[c] = max(maxima[c], best)
        counts[c] += 1
    return maxima, counts


def recurrence_gaps(a, N):
    u, v, out = 0, 1, []
    while v < N:
        out.append(v)
        u, v = v, a*v+u
    return out


def main():
    witness = proper(8, 4)
    assert witness is not None and proper(9, 4) is None
    assert all(witness[y-1] != witness[x-1]
               for x in range(1, 9) for y in range(x+1, 9)
               if y-x in fib_gaps(8))
    print('FIDELITY Delta(F,2;4)=9; avoiding colors on [8]:', witness)
    cases = 0
    for a in range(4, 41, 4):
        D, lam, alpha = a*a+4, (Q(a, 2), Q(-1, 2)), (Q(1-a, 4), Q(1, 4))
        assert sign(lam, D) < 0 and sign(add(lam, (Q(1, 4), Q(0))), D) > 0
        assert mul(lam, lam, D) == add(scale(a, lam), (Q(1), Q(0)))
        u, v, h, j, power = 0, 1, 2, 1, (Q(1), Q(0))
        for n in range(65):
            assert add(scale(4*u, alpha), scale(2, power)) == (Q(h), Q(0))
            assert h % 4 == (1 if n % 2 else 2)
            if n:
                phase = add(scale(u, alpha), (-Q(floor_alpha(a, u)), Q(0)))
                expected = add((Q(1, 4) if n % 2 else Q(1, 2), Q(0)), scale(Q(-1, 2), power))
                assert phase == expected
                assert sign(add(phase, (Q(-1, 4), Q(0))), D) > 0
                assert sign(add(phase, (Q(-1, 2), Q(0))), D) < 0
                cases += 1
            power = mul(power, lam, D)
            u, v, h, j = v, a*v+u, j, a*j+h
    print('EXACT_ALGEBRA a=4,8,...,40; n=1,...,64; identities and phase bounds:', cases)
    f = [0, 1]
    for _ in range(605):
        f.append(f[-1]+f[-2])
    u, v = 0, 1
    for n in range(201):
        assert 2*u == f[3*n]
        assert (f[n] % 2 == 0) == (n % 3 == 0)
        assert f[n+6] == 4*f[n+3]+f[n]
        u, v = v, 4*v+u
    print('FIBONACCI u_n=F_(3n)/2, parity, six-step recurrence: n=0,...,200')
    gaps = recurrence_gaps(4, 50000)
    maxima, counts = longest(50000, gaps, lambda n: chi(4, n), 2)
    assert maxima == [2, 2]
    print('FIDELITY Delta(G,3;2)>50000; all gaps:', gaps)
    print('G colors on [50000]; counts:', counts, 'longest lengths:', maxima)
    gaps = fib_gaps(99999)
    maxima, counts = longest(100000, gaps, lambda n: 2*(n % 2)+chi(4, n//2+1), 4)
    assert maxima == [2, 2, 2, 2]
    print('SETTLEMENT_CHECK F on [100000]; distinct Fibonacci gaps:', len(gaps), gaps)
    print('F color counts:', counts, 'longest lengths:', maxima)
    for a in (4, 8, 12, 16):
        gaps = recurrence_gaps(a, 10000)
        maxima, _ = longest(10000, gaps, lambda n, a=a: chi(a, n), 2)
        assert max(maxima) <= 2
        print('FAMILY_CHECK a=', a, 'N=10000; gaps:', gaps, 'longest lengths:', maxima)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
