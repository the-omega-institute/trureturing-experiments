#!/usr/bin/env python3
"""Exact checks for Hardin's OEIS A250538--A250544 empirical recurrences.

NumPy is used only as a container for arbitrary-precision Python integers
(dtype=object). No floating-point arithmetic or fitted recurrence is used.
"""
from collections import Counter
from itertools import product
from pathlib import Path
import json
import re
import numpy as np

HERE = Path(__file__).resolve().parent
SOURCE = json.loads((HERE / 'source-data.json').read_text())
COEFF = [1, -30, 415, -3514, 20386, -85924, 272198, -661180,
         1244717, -1822478, 2068955, -1802474, 1181760, -563888,
         184752, -37152, 3456]
SHORT = [1, -16, 106, -376, 769, -904, 564, -144]


def bfile(a):
    return {int(n): int(v) for line in SOURCE[a]['bfile'].splitlines()
            if line.strip() and not line.startswith('#')
            for n, v in [line.split()]}


def polynomial(factors):
    c = [1]
    for w, m in factors:
        for _ in range(m):
            d = [0] * (len(c) + 1)
            for i, v in enumerate(c):
                d[i] += v
                d[i + 1] -= w * v
            c = d
    return c


def row_states(width):
    # Independent enumeration of the ORIGINAL possible rows.
    c = Counter(tuple(r[j + 1] - r[j] for j in range(width - 1))
                for r in product(range(4), repeat=width))
    for d, w in c.items():
        p = [0]
        for v in d:
            p.append(p[-1] + v)
        assert w == 4 - max(p) + min(p)
        assert -3 <= sum(d) <= 3
    return c


def counts(width, max_n):
    states = row_states(width)
    dim = width - 1
    keys = list(states)
    inds = tuple(np.array([d[j] + 3 for d in keys], dtype=int)
                 for j in range(dim))
    weights = np.array([states[d] for d in keys], dtype=object)
    vals = weights.copy()
    out = {0: sum(vals)}  # one row, the natural n=0 extension
    for n in range(1, max_n + 1):
        z = np.zeros((7,) * dim, dtype=object)
        z[inds] = vals
        # Product-poset zeta transform: sum all d <= e, including equality.
        for axis in range(dim):
            z = z.cumsum(axis=axis, dtype=object)
        vals = weights * z[inds]
        out[n] = sum(vals)
        assert type(out[n]) is int
    print(f'COUNT width={width} states={len(states)} n=0..{max_n} '
          f'first={out[1]} last={out[max_n]}', flush=True)
    return out


def direct(width, max_n):
    # Actual rows with all 4^width offsets retained, and the source's two
    # monotonicity conditions checked separately on consecutive rows.
    rows = list(product(range(4), repeat=width))
    def valid(x, y):
        horizontal = all(x[j + 1] - x[j] <= y[j + 1] - y[j]
                         for j in range(width - 1))
        vertical = all(y[j] - x[j] <= y[j + 1] - x[j + 1]
                       for j in range(width - 1))
        assert horizontal == vertical
        return horizontal and vertical
    incoming = [[i for i, x in enumerate(rows) if valid(x, y)] for y in rows]
    v = [1] * len(rows)
    result = {0: len(rows)}
    for n in range(1, max_n + 1):
        v = [sum(v[i] for i in inc) for inc in incoming]
        result[n] = sum(v)
    return result


def residuals(values, coeff, start):
    for n in range(start, max(values) + 1):
        assert sum(c * values[n - j] for j, c in enumerate(coeff)) == 0, n


def matrix_check(width, coeff):
    states = row_states(width)
    ds = sorted(states)
    T = [[states[e] if all(x <= y for x, y in zip(d, e)) else 0
          for e in ds] for d in ds]
    m = len(ds)
    # Horner evaluation of t^r + c1*t^(r-1) + ... + cr at T.
    M = [[int(i == j) for j in range(m)] for i in range(m)]
    for c in coeff[1:]:
        M = [[sum(M[i][k] * T[k][j] for k in range(m)) +
              (c if i == j else 0) for j in range(m)] for i in range(m)]
    assert all(v == 0 for row in M for v in row)
    print(f'MATRIX width={width} size={m} polynomial_degree={len(coeff)-1} zero=True')


def main():
    assert polynomial([(1, 7), (2, 5), (3, 3), (4, 1)]) == COEFF
    assert polynomial([(1, 2), (2, 2), (3, 2), (4, 1)]) == SHORT
    print('DENOMINATORS exact coefficient vectors match both printed recurrences')
    computed = {2: counts(2, 32)}
    for a in range(250538, 250544):
        name = f'A{a}'
        source = SOURCE[name]['text']
        formula = next(x for x in source.splitlines() if x.startswith('%F'))
        terms = [(int(c.replace(' ', '')), int(j)) for c, j in
                 re.findall(r'([+-]?\s*\d+)\*a\(n-(\d+)\)', formula)]
        assert terms == [(-COEFF[j], j) for j in range(1, 17)]
        assert f'%O {name} 1,1' in source and 'Empirical:' in formula
        width = a - 250535
        expected = bfile(name)
        computed[width] = counts(width, max(expected))
        assert all(computed[width][n] == v for n, v in expected.items())
        displayed = [int(v) for line in source.splitlines()
                     if re.match(r'%[STU] ', line)
                     for v in line.split(' ', 2)[2].split(',') if v]
        assert displayed == [expected[n] for n in range(1, len(displayed) + 1)]
        residuals(computed[width], COEFF, 16)
        print(f'FIDELITY {name} bfile_terms={len(expected)} '
              f'n={min(expected)}..{max(expected)} displayed={len(displayed)} '
              f'recurrence_n=16..{max(expected)} OK', flush=True)
    for width in [2, 3, 4]:
        d = direct(width, 24)
        assert all(v == computed[width][n] for n, v in d.items())
        print(f'INDEPENDENT actual-row DP width={width} n=0..24 agrees')
    source = SOURCE['A250544']['text']
    f = next(x for x in source.splitlines() if x.startswith('%F') and 'k=1:' in x)
    terms = [(int(c.replace(' ', '')), int(j)) for c, j in
             re.findall(r'([+-]?\s*\d+)\*a\(n-(\d+)\)', f)]
    assert terms == [(-SHORT[j], j) for j in range(1, 8)]
    residuals(computed[2], SHORT, 7)
    residuals(computed[2], COEFF, 16)
    table = bfile('A250544')
    checked = 0
    n, k = 1, 1
    for index, value in sorted(table.items()):
        assert index >= 1
        if k + 1 in computed:
            assert computed[k + 1][n] == value
            checked += 1
        elif n + 1 in computed:
            # Transposition preserves BOTH literal source conditions.
            assert computed[n + 1][k] == value
            checked += 1
        if k == 1:
            k, n = n + 1, 1
        else:
            k, n = k - 1, n + 1
    assert checked > 130
    print(f'FIDELITY A250544 antidiagonal bfile matched={checked}/{len(table)}; '
          'remaining entries outside width<=8 are explicitly untested')
    print('A250544 column1 order7 n=7..32 and order16 n=16..32 OK; '
          'columns2..7 checked above')
    matrix_check(2, SHORT)
    matrix_check(3, COEFF)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
