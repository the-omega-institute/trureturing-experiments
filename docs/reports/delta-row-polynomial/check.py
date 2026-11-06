#!/usr/bin/env python3
"""Exact checks for OEIS A111184, Kurkov Conjecture 3.

Polynomials are tuples of integer coefficients in ascending degree.
The DELTA implementation retains both variables before evaluating X=1.
These finite checks are evidence for fidelity and indexing, not a proof.
"""
from functools import lru_cache
from math import factorial


def norm(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return tuple(a)


def add(a, b):
    c = [0] * max(len(a), len(b))
    for i, v in enumerate(a):
        c[i] += v
    for i, v in enumerate(b):
        c[i] += v
    return norm(c)


def neg(a):
    return tuple(-x for x in a)


def mul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return norm(c)


def evaluate(a, q):
    out = 0
    for x in reversed(a):
        out = q * out + x
    return out


@lru_cache(None)
def delta(n, k):
    """Literal bivariate p(n,k), represented by {(X-degree,Y-degree): integer}."""
    if n == 0:
        return {(0, 0): 1}
    if k == -1:
        return {}
    r = k // 2 if k % 2 == 0 else k // 2 + 2
    s = 1 if k % 2 == 0 else 0
    out = dict(delta(n, k - 1))
    for (i, j), v in delta(n - 1, k + 1).items():
        out[i + 1, j] = out.get((i + 1, j), 0) + r * v
        out[i, j + 1] = out.get((i, j + 1), 0) + s * v
    return {ij: v for ij, v in out.items() if v}


def row(n):
    out = [0] * (n + 1)
    for (i, j), v in delta(n, 0).items():
        assert i + j == n
        out[j] += v
    return norm(out)


def loop(m):
    v = [(1,)] * (m + 1)  # entry zero is unused
    for i in range(1, m):
        acc = v[i]
        for j in range(i + 1, m + 1):
            old = v[j]
            acc = add(acc, old)
            v[j] = add(mul((j - i - 1, 1), old), acc)
    return v


def quotient(N):
    h, hp = [(1,)], [(1,)]
    for n in range(1, N + 1):
        h.append(mul(h[-1], (n - 1, 1)))
        hp.append(mul(hp[-1], (n, 1)))
    out = []
    for n in range(N + 1):
        c = hp[n]
        for j in range(1, n + 1):
            c = add(c, neg(mul(h[j], out[n - j])))
        out.append(c)
    return out


def return_walks(N):
    state, out = {0: (1,)}, [(1,)]
    for _ in range(N):
        nxt = {}
        for k, v in state.items():
            for s in range(k + 1):
                nxt[s] = add(nxt.get(s, (0,)), v)
            nxt[k + 1] = add(nxt.get(k + 1, (0,)), mul((k + 1, 1), v))
        state = nxt
        out.append(state[0])
    return out


def main():
    # OEIS A111184 %e, internal revision 58: all displayed complete rows.
    published = [
        [1], [0, 1], [0, 2, 1], [0, 6, 6, 1],
        [0, 24, 34, 12, 1], [0, 120, 210, 110, 20, 1],
        [0, 720, 1452, 974, 270, 30, 1],
        [0, 5040, 11256, 8946, 3248, 560, 42, 1],
        [0, 40320, 97296, 87504, 38338, 8792, 1036, 56, 1],
    ]
    for n, expected in enumerate(published):
        assert row(n) == tuple(expected), (n, row(n), expected)
    print('FIDELITY: all 9 displayed A111184 rows n=0..8 reproduced exactly')
    N = 20
    series, walks = quotient(N), return_walks(N)
    count = 0
    for m in range(1, N + 1):
        v = loop(m)
        for n in range(1, m + 1):
            assert row(n) == mul((0, 1), v[n]), (m, n, 'DELTA versus loop')
            assert v[n] == series[n] == walks[n], (m, n, 'quotient / walks')
            count += 1
    print(f'IDENTITY: all {count} pairs 1 <= n <= m <= {N}, in Z[q]')
    print(f'INDEPENDENT BRIDGES: quotient H_(q+1)/H_q and weighted return walks through n={N}')
    for n in range(1, 6):
        print(f'v_{n}, coefficients ascending: {series[n]}')
    for n in range(1, N + 1):
        assert evaluate(series[n], 0) == factorial(n)
        for q in range(-5, 7):
            assert evaluate(row(n), q) == q * evaluate(series[n], q)
    print('SPECIALIZATIONS: q=-5..6, n=1..20; at q=0, v_n=n!')
    # An independent scalar U recurrence, already owned by the repository.
    @lru_cache(None)
    def U(n, k):
        return 1 if n == 0 else (k + 2) * U(n - 1, k + 1) + sum(U(n - 1, j) for j in range(k + 1))
    for n in range(1, N + 1):
        assert evaluate(series[n], 1) == U(n - 1, 0)
    print('REPOSITORY SCALAR BRIDGE: v_n(1)=U(n-1,0), n=1..20')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
