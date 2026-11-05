#!/usr/bin/env python3
"""Exact finite checks for Hanna's OEIS A193550 square-root conjecture.

Python 3 standard library only. Coefficients are in ascending degree order.
The long computation uses weighted Dyck paths; an independent direct nested
fraction computation checks the first 41 coefficients. Neither is an all-order
proof. Source: https://oeis.org/A193550/internal, revision 19.
"""
from fractions import Fraction
from math import comb


SOURCE_A = [
    1, 1, 2, 13, 206, 5794, 252068, 15663997, 1316280854,
    143694972886, 19764465724412, 3343418236081618,
    682133942067492236, 165161123584687819684, 46817735849074712020808,
    15358634840651231221695517, 5772973821383087169122348774,
]
SOURCE_B = [
    1, 1, 6, 97, 2782, 122670, 7687932, 649446621,
    71136143478, 9806113041658,
]


def beta(j):
    return (2 * (j // 2) + 1) ** 2


def mul(f, g, n):
    return [sum(f[i] * g[k - i] for i in range(k + 1))
            for k in range(n + 1)]


def inverse(f):
    assert f[0] == 1
    g = [1]
    for k in range(1, len(f)):
        g.append(-sum(f[i] * g[k - i] for i in range(1, k + 1)))
    return g


def fraction(n, depth, weights=beta, terminal=None):
    """T_0(...T_(depth-1)(terminal)...), truncated after degree n."""
    f = [1] + [0] * n if terminal is None else terminal[:]
    assert len(f) == n + 1 and f[0] == 1
    for j in reversed(range(depth)):
        f = inverse([1] + [-weights(j) * f[k - 1]
                           for k in range(1, n + 1)])
    return f


def dyck_moments(n, weights=beta):
    """Sum paths staying nonnegative; down h -> h-1 has weight beta(h-1)."""
    dp = [0] * (n + 1)
    dp[0] = 1
    a = [1]
    for t in range(1, 2 * n + 1):
        nxt = [0] * (n + 1)
        for h, v in enumerate(dp):
            if h < n:
                nxt[h + 1] += v
            if h:
                nxt[h - 1] += weights(h - 1) * v
        dp = nxt
        if t % 2 == 0:
            a.append(dp[0])
    return a


def normalized_integral_root(f):
    assert f[0] == 1
    b = [1]
    for k in range(1, len(f)):
        numerator = f[k] - sum(b[i] * b[k - i] for i in range(1, k))
        q, remainder = divmod(numerator, 2)
        assert remainder == 0, ("nonintegral root", k, Fraction(numerator, 2))
        b.append(q)
    assert mul(b, b, len(f) - 1) == f
    return b


def compose(f, g, n):
    assert g[0] == 0
    out = [0] * (n + 1)
    for coefficient in reversed(f[:n + 1]):
        out = mul(out, g, n)
        out[0] += coefficient
    return out


def main():
    n = 500
    a = dyck_moments(n)
    assert a[:len(SOURCE_A)] == SOURCE_A
    f = a[1:]
    b = normalized_integral_root(f)
    assert b[:len(SOURCE_B)] == SOURCE_B
    print("FIDELITY: all 17 displayed A coefficients and 10 displayed B coefficients match")
    print("FIDELITY: B coefficients 0..499 are integers, reproducing the source's initial 500 terms")
    print("B[0..9] =", b[:10])

    m = 40
    direct = fraction(m, m)
    assert direct == a[:m + 1]
    assert direct == fraction(m, m + 5)
    terminal = [1, 3, 7] + [0] * (m - 2)
    assert direct == fraction(m, m, terminal=terminal)
    assert dyck_moments(m, lambda j: 1) == [comb(2*k, k)//(k+1) for k in range(m+1)]
    print("INDEPENDENCE: direct finite fraction matches weighted-path computation through degree 40")
    print("STABILIZATION: depths 40 and 45 and a nontrivial normalized terminal tail agree through degree 40")

    c = [comb(2*k, k)//(k+1) for k in range(n + 1)]
    assert all((a[k] - c[k]) % 8 == 0 for k in range(n + 1))
    cc = mul(c[:n], c[:n], n - 1)
    assert c[1:] == cc
    assert all((f[k] - cc[k]) % 8 == 0 for k in range(n))
    print("MODULAR: A = Catalan modulo 8 through degree 500; F = Catalan^2 modulo 8 through degree 499")

    # Test the proposed explicit lift independently of the root recurrence.
    # All operations below are exact truncated integer-series operations.
    h = [(f[k] - cc[k]) // 8 for k in range(m + 1)]
    assert h[0] == 0
    cinv = inverse(c[:m + 1])
    u = [2*v for v in mul(h, mul(cinv, cinv, m), m)]
    d = compose(c[:m + 1], [-v for v in u], m)
    ud = mul(u, d, m)
    ud2 = mul(u, mul(d, d, m), m)
    assert [d[k] + ud2[k] for k in range(m + 1)] == [1] + [0]*m
    v = [1 + 2*ud[0]] + [2*q for q in ud[1:]]
    assert mul(v, v, m) == [1 + 4*u[0]] + [4*q for q in u[1:]]
    lifted = mul(c[:m + 1], v, m)
    assert lifted == b[:m + 1]
    assert mul(lifted, lifted, m) == f[:m + 1]
    print("LIFT: D + U*D^2 = 1, V^2 = 1 + 4*U, and (C*V)^2 = F through degree 40")
    print("EXPECTED SETTLEMENT CHECK: no failure in the tested finite ranges; no all-order claim is certified")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
