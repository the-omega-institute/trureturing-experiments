"""Exact finite checks for the tropical TN coefficient-space closure proposal.

Python 3 standard library only. None denotes minus infinity. These checks
test the general argument; no finite enumeration proves the universal claim.
"""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import factorial


def sign(p):
    return (-1) ** sum(p[i] > p[j] for i in range(len(p))
                       for j in range(i + 1, len(p)))


def weight(a, rows, cols, p):
    entries = [a[i][cols[p[r]]] for r, i in enumerate(rows)]
    return None if None in entries else sum(entries)


def leq(x, y):
    return x is None or (y is not None and x <= y)


def plus(x, y):
    return None if x is None or y is None else x + y


def monge(a):
    n = len(a)
    return all(leq(plus(a[i][l], a[j][k]), plus(a[i][k], a[j][l]))
               for i, j in combinations(range(n), 2)
               for k, l in combinations(range(n), 2))


def minors(n):
    return [(r, c, tuple(permutations(range(k))))
            for k in range(1, n + 1)
            for r in combinations(range(n), k)
            for c in combinations(range(n), k)]


def dominant(a, rows, cols, ps):
    ws = [weight(a, rows, cols, p) for p in ps]
    finite = [w for w in ws if w is not None]
    if not finite:
        return None, []
    w = max(finite)
    return w, [p for p, wp in zip(ps, ws) if wp == w]


def delta(c, rows, cols, opts):
    out = 0
    for p in opts:
        term = sign(p)
        for r, i in enumerate(rows):
            term *= c[i][cols[p[r]]]
        out += term
    return out


def check_matrix(a, data):
    n = len(a)
    b = 1 + 3 ** n * factorial(n)
    c = [[0 if a[i][j] is None else b ** ((i + 1) * (j + 1))
          for j in range(n)] for i in range(n)]
    checked = tied = absent = 0
    for rows, cols, ps in data:
        w, opts = dominant(a, rows, cols, ps)
        if w is None:
            absent += 1
            continue
        k = len(rows)
        ordered = tuple(range(k))
        assert ordered in opts, ('uncrossing', a, rows, cols)
        sid = sum((i + 1) * (j + 1) for i, j in zip(rows, cols))
        for p in opts:
            if p != ordered:
                sp = sum((rows[r] + 1) * (cols[p[r]] + 1) for r in range(k))
                assert sp <= sid - 1
        # This exact interval bound holds for EVERY point in the open box.
        negatives = sum(sign(p) < 0 for p in opts)
        lower = Q(b ** sid, 2 ** k)
        upper = negatives * Q(3 ** k * b ** (sid - 1), 2 ** k)
        assert lower > upper
        assert b > 3 ** k * (factorial(k) - 1)
        assert delta(c, rows, cols, opts) > 0
        checked += 1
        tied += len(opts) > 1
    return checked, tied, absent


def exhaustive(n, diagonal_values, off_values):
    data = minors(n)
    off = [(i, j) for i in range(n) for j in range(n) if i != j]
    total = accepted = checked = tied = absent = 0
    for ds in product(diagonal_values, repeat=n):
        for es in product(off_values, repeat=len(off)):
            a = [[None] * n for _ in range(n)]
            for i in range(n):
                a[i][i] = ds[i]
            for (i, j), value in zip(off, es):
                a[i][j] = value
            total += 1
            if not monge(a):
                continue
            accepted += 1
            x, y, z = check_matrix(a, data)
            checked += x
            tied += y
            absent += z
    print(f'EXHAUSTIVE n={n} total={total} tropical_TN={accepted} '
          f'finite_minors={checked} tied_minors={tied} zero_minors={absent}')
    return accepted, checked


def padd(a, b):
    out = dict(a)
    for e, c in b.items():
        out[e] = out.get(e, 0) + c
        if out[e] == 0:
            del out[e]
    return out


def pmul(a, b):
    out = {}
    for e, c in a.items():
        for f, d in b.items():
            out[e + f] = out.get(e + f, 0) + c * d
    return {e: c for e, c in out.items() if c != 0}


def det(a):
    out = {}
    for p in permutations(range(len(a))):
        term = {Q(0): sign(p)}
        for i in range(len(a)):
            term = pmul(term, a[i][p[i]])
        out = padd(out, term)
    return out


def term(e, c=1):
    return {Q(e): Q(c)} if c else {}


def fidelity():
    a = [[10, 9, 4], [9, 8, 4], [4, 4, 2]]
    assert monge(a)
    for r, s in [(Q(2), Q(1)), (Q(5, 3), Q(2, 3)), (Q(7, 4), Q(7, 4))]:
        m = [[term(10), term(9, s), term(4)],
             [term(9), term(8, r), term(4)],
             [term(4), term(4), term(2)]]
        expected = padd(padd(term(20, r-s), term(18, -1)),
                        padd(term(17, s+1), term(16, -r)))
        assert det(m) == expected
        c = [[1, s, 1], [1, r, 1], [1, 1, 1]]
        ds = []
        for k in (1, 2, 3):
            terms = []
            for rows in combinations(range(3), k):
                w, opts = dominant(a, rows, rows, tuple(permutations(range(k))))
                terms.append((w, delta(c, rows, rows, opts)))
            maximum = max(w for w, _ in terms)
            ds.append(sum(d for w, d in terms if w == maximum))
        assert ds == [1, r-s, r-s]
    for eps in (Q(1, 2), Q(1), Q(3, 2)):
        m = [[term(10), term(9), term(4)],
             [term(9), padd(term(8), term(8-eps)), term(4)],
             [term(4), term(4), term(2)]]
        expected = padd(padd(term(20-eps), term(18, -1)),
                        padd(term(17, 2), padd(term(16, -1), term(16-eps, -1))))
        assert det(m) == expected
        assert max(det(m)) == 20-eps
        for rows, cols, _ in minors(3):
            d = det([[m[i][j] for j in cols] for i in rows])
            assert not d or d[max(d)] > 0
    b = [[1, 2, 3], [2, 1, 4], [3, 4, 1]]
    c = [[2, 3, 1], [4, 5, 6], [7, 8, 9]]
    w, opts = dominant(b, tuple(range(3)), tuple(range(3)), tuple(permutations(range(3))))
    assert w == 9 and delta(c, range(3), range(3), opts) == 62
    print('FIDELITY Example 5.3: permanent=9 dominant_coefficient=62')
    print('FIDELITY Example 5.17: D=(1,r-s,r-s), exact determinants, epsilon cancellations OK')


def tied_first():
    a = [[0, 0], [0, 0]]
    _, opts = dominant(a, (0, 1), (0, 1), tuple(permutations(range(2))))
    assert delta([[1, 1], [1, 1]], (0, 1), (0, 1), opts) == 0
    assert delta([[1, 1], [1, 2]], (0, 1), (0, 1), opts) == 1
    assert det([[term(0), term(0)], [padd(term(0), term(-1,-1)),
                padd(term(0), term(-1))]]) == term(-1, 2)
    print('TIE_FIRST: same valuation A=0 admits Delta=0 and Delta=1; tie does not force equality')


def box_vertices():
    # Exact multiaffine corner evaluations certify these entire closed boxes.
    a = [[0] * 3 for _ in range(3)]
    b = 1 + 3 ** 3 * factorial(3)
    minimum = None
    for choices in product((Q(1, 2), Q(3, 2)), repeat=9):
        c = [[choices[3*i+j] * b ** ((i+1)*(j+1)) for j in range(3)] for i in range(3)]
        for rows, cols, ps in minors(3):
            _, opts = dominant(a, rows, cols, ps)
            d = delta(c, rows, cols, opts)
            assert d > 0
            minimum = d if minimum is None else min(minimum, d)
    print(f'BOX_VERTICES n=3 vertices=512 all_19_minors_positive minimum={minimum}')


def main():
    tied_first()
    fidelity()
    for n in (1, 2, 3):
        exhaustive(n, (-1, 0, 1), (None, -1, 0, 1))
    exhaustive(4, (0,), (None, 0))
    box_vertices()
    print('SCOPE finite exact checks and full-box interval bounds; universal theorem remains unproved in Lean')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
