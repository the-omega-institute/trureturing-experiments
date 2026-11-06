#!/usr/bin/env python3
"""Exact, independent checks for Batavia--Javadekar--Kumar Q6.2.

Only Python's integer arithmetic is used. Finite-field elimination is explicit.
The source is arXiv:2609.26752v1, Theorem A(b,c), Question 6.2.
This program checks finite instances; it does not prove the infinite claim.
"""

from itertools import combinations, product
from math import comb


def choose(n, r):
    return comb(n, r) if 0 <= r <= n else 0


def rank(columns, p):
    """Rank of sparse columns over the prime field F_p."""
    pivots = {}
    for column in columns:
        v = {i: x % p for i, x in column.items() if x % p}
        while v:
            i = min(v)
            if i not in pivots:
                a = pow(v[i], -1, p)
                pivots[i] = {j: a * x % p for j, x in v.items()}
                break
            a = v[i]
            for j, x in pivots[i].items():
                y = (v.get(j, 0) - a * x) % p
                if y:
                    v[j] = y
                else:
                    v.pop(j, None)
    return len(pivots)


def apery(m, n):
    """Compute Ap(S,m) from the generating set, without the proposed table."""
    gens = [g for g in range(m, 2 * m) if g not in (m + 2, m + n)]
    reachable = {0}
    a = [None] * m
    for x in range(5 * m):
        if x and any(x - g in reachable for g in gens):
            reachable.add(x)
        if x in reachable and a[x % m] is None:
            a[x % m] = x
    assert None not in a
    expected = [0] + [2 * m + r if r in (2, n) else m + r
                      for r in range(1, m)]
    assert a == expected, (m, n, a, expected)
    positive = a[1:]
    index = {x: i for i, x in enumerate(positive)}
    table = {(i, j): index.get(x + y)
             for i, x in enumerate(positive) for j, y in enumerate(positive)}
    # Independent comparison with the claimed residue-indexed multiplication.
    allowed = set(range(1, m)) - {2, n}
    for i in range(1, m):
        for j in range(1, m):
            want = None
            if i in allowed and j in allowed:
                if i == j == 1:
                    want = 1  # residue 2 has index 1 in positive
                elif i + j == n:
                    want = n - 1
            assert table[i - 1, j - 1] == want, (m, n, i, j)
    return a, table


def add(v, i, a):
    v[i] = v.get(i, 0) + a
    if not v[i]:
        del v[i]


def bar_column(word, table):
    v = {}
    for i in range(len(word) - 1):
        t = table[word[i], word[i + 1]]
        if t is not None:
            w = word[:i] + (t,) + word[i + 2:]
            add(v, w, (-1) ** i)
    return v


def bar_betti(m, n, p, degree=4):
    """Reduced normalized bar homology of A over F_p through degree 4."""
    a, table = apery(m, n)
    r = m - 1
    ranks = [0, 0]
    for d in range(2, degree + 2):
        def columns():
            for word in product(range(r), repeat=d):
                v = bar_column(word, table)
                # Check d^2=0 integrally, independently of the rank computation.
                vv = {}
                for w, c in v.items():
                    for u, b in bar_column(w, table).items():
                        add(vv, u, c * b)
                assert not vv
                yield v
        ranks.append(rank(columns(), p))
    actual_a = [1] + [r ** d - ranks[d] - ranks[d + 1]
                      for d in range(1, degree + 1)]
    e = m - 3
    delta = int(n >= 6)
    expected_a = [1, e]
    for i in range(2, degree + 1):
        expected_a.append(e * expected_a[-1] - delta * expected_a[-2])
    assert actual_a == expected_a, (m, n, p, actual_a, expected_a)
    actual_r = [1] + [actual_a[d] + actual_a[d - 1]
                      for d in range(1, degree + 1)]
    print(f"BAR m={m} n={n} F_{p}: ranks(d_1..d_{degree+1})={ranks[1:]} "
          f"beta_A={actual_a}; regular-element transfer predicts beta_R={actual_r}")
    return actual_r


def koszul_betti(m, n, p):
    """Koszul homology on A's degree-one generators; Tor^Q(R,k)."""
    a, table = apery(m, n)
    residues = sorted(set(range(1, m)) - {2, n})
    e = len(residues)
    ranks = [0]
    for d in range(1, e + 1):
        def columns():
            for subset in combinations(range(e), d):
                for b in range(m):  # 0 is 1; i>0 is residue i
                    v = {}
                    for j, g in enumerate(subset):
                        residue = residues[g]
                        t = residue if b == 0 else table[b - 1, residue - 1]
                        if t is not None:
                            if b:
                                t += 1
                            add(v, (t, subset[:j] + subset[j+1:]), (-1) ** j)
                    yield v
        ranks.append(rank(columns(), p))
    ranks.append(0)
    actual = [m * choose(e, d) - ranks[d] - ranks[d + 1]
              for d in range(e + 1)]
    expected = [1] + [d * choose(m - 2, d + 1)
        - (choose(m - 3, d - 1) if n >= 6 else 0)
        + (choose(m - n + 2, d - n + 4) if n >= 6 else 0)
        for d in range(1, e + 1)]
    assert actual == expected, (m, n, p, actual, expected)
    print(f"FIDELITY Theorem A({'c' if n >= 6 else 'b'}), "
          f"m={m} n={n} F_{p}: B_coeff={actual}")
    return actual


def multiply(a, b, limit):
    out = [0] * (limit + 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j <= limit:
                out[i + j] += x * y
    return out


def inverse(a, limit):
    assert a[0] == 1
    out = [1]
    for i in range(1, limit + 1):
        out.append(-sum(a[j] * out[i - j]
                        for j in range(1, min(i, len(a) - 1) + 1)))
    return out


def series_checks(m, n, limit=30):
    e = m - 3
    delta = int(n >= 6)
    b = [1] + [d * choose(m - 2, d + 1)
        - (choose(m - 3, d - 1) if delta else 0)
        + (choose(m - n + 2, d - n + 4) if delta else 0)
        for d in range(1, e + 1)]
    denominator = [1, 0] + [-x for x in b[1:]]
    serre = multiply([choose(e + 1, i) for i in range(e + 2)],
                     inverse(denominator, limit), limit)
    q = [1, -e, delta]
    poincare = multiply([1, 1], inverse(q, limit), limit)
    defect = [x - y for x, y in zip(serre, poincare)]
    assert all(x >= 0 for x in defect)
    if not delta:
        assert not any(defect)
    else:
        assert next(i for i, x in enumerate(defect) if x) == n - 3
        assert defect[n - 3] == 1
        # Check the entire asserted difference identity, not only its order.
        t = multiply([choose(n - 5, i) for i in range(n - 4)], q, limit)
        t[n - 3] -= 1
        rhs = multiply(inverse(q, limit), inverse(t, limit), limit)
        rhs = multiply([0] * (n - 3) + [1, 1], rhs, limit)
        assert defect == rhs, (m, n)


def dual_dimension(m, n, p, d):
    """Degree-d quotient rank using only the TWO defining dual relations."""
    alphabet = sorted(set(range(1, m)) - {2, n})
    rels = [{(1, 1): 1}, {(i, n - i): 1 for i in alphabet
                                        if n - i in alphabet}]
    def columns():
        for r in rels:
            for left_len in range(d - 1):
                for left in product(alphabet, repeat=left_len):
                    for right in product(alphabet, repeat=d - 2 - left_len):
                        yield {left + w + right: c for w, c in r.items()}
    return len(alphabet) ** d - rank(columns(), p)


def normal_count(m, n, d):
    """Finite-state enumeration of words avoiding the proposed leading words."""
    alphabet = sorted(set(range(1, m)) - {2, n})
    q = [(i, n - i) for i in alphabet if n - i in alphabet
         and i not in (1, n - 1)]
    # a>b>others; among the others natural order determines degree-lex.
    order = {a: (2 if a == 1 else 1 if a == n - 1 else 0, a)
             for a in alphabet}
    leading = max(q, key=lambda w: tuple(order[x] for x in w)) if q else None
    forbidden = [(1, 1), (1, n - 1)] + ([(1,) + leading] if leading else [])
    states = {(): 1}
    for _ in range(d):
        new = {}
        for suffix, c in states.items():
            for a in alphabet:
                word = suffix + (a,)
                if any(len(word) >= len(w) and word[-len(w):] == w for w in forbidden):
                    continue
                new[word[-2:]] = new.get(word[-2:], 0) + c
        states = new
    return sum(states.values())


def main():
    print("Exact arithmetic: integers and modular inverses in F_2, F_3; no floating point")
    for m, n in [(5, 4), (6, 5), (7, 6)]:
        for p in (2, 3):
            koszul_betti(m, n, p)
            bar_betti(m, n, p)
    count = 0
    for m in range(5, 16):
        for n in range(4, m):
            apery(m, n)
            series_checks(m, n)
            count += 1
    print(f"APERY_AND_SERRE: {count} pairs, 5<=m<=15, 4<=n<m; series through degree 30")
    for m, n in [(5, 4), (6, 5), (7, 6), (8, 7), (9, 8)]:
        e = m - 3
        # (1+z)/H = 1+(1-e)z+(2-e)z^2+delta*z^3
        den = [1, 1-e, 2-e, int(n >= 6)]
        h = multiply([1, 1], inverse(den, 6), 6)
        for d in range(7):
            assert normal_count(m, n, d) == h[d]
        for p in (2, 3):
            for d in range(2, 5):
                actual = dual_dimension(m, n, p, d)
                assert actual == h[d], (m, n, p, d, actual, h[d])
        print(f"DUAL m={m} n={n}: H[0..6]={h}; two-relation rank verified d=2..4 in F_2,F_3")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
