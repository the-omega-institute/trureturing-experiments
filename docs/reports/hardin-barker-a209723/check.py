#!/usr/bin/env python3
"""Exact, finite fidelity and conjecture checks for OEIS A209723.

Source data: OEIS Foundation and contributors, b209723.txt, CC BY-SA 4.0.
These finite checks are not a proof of an assertion at unbounded height.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
from pathlib import Path


def square_ok(a, b, c, d):
    return len({b - a, d - b, c - d, a - c}) == 4


def compatible(s, t):
    return all(square_ok(s[j], s[j + 1], t[j], t[j + 1]) for j in range(4))


def digits(s):
    return tuple(map(int, s))


def direct_count(n):
    """Place individual entries, pruning only a completed literal square.

    This does not use the row graph, its component list, or the parity formulas.
    """
    cells = []

    def visit():
        k = len(cells)
        if k == 5 * (n + 1):
            return 1
        r, c = divmod(k, 5)
        result = 0
        for v in range(3):
            cells.append(v)
            valid = True
            if r > 0 and c > 0:
                cycle = (cells[k - 6], cells[k - 5], cells[k], cells[k - 1])
                delta = [cycle[(i + 1) % 4] - cycle[i] for i in range(4)]
                valid = all(delta[i] != delta[j] for i in range(4) for j in range(i))
            if valid:
                result += visit()
            cells.pop()
        return result

    return visit()


def main():
    squares = {v for v in product(range(3), repeat=4) if square_ok(*v)}
    expected_squares = {
        (0, 1, 2, 0), (0, 2, 1, 0), (0, 2, 2, 1), (1, 0, 0, 2),
        (1, 2, 2, 0), (2, 0, 0, 1), (2, 0, 1, 2), (2, 1, 0, 2),
    }
    assert squares == expected_squares
    column_arcs = {((a, c), (b, d)) for a, b, c, d in squares}
    expected_arcs = {(digits(a), digits(b)) for a, b in [
        ('02', '10'), ('02', '21'), ('20', '01'), ('20', '12'),
        ('10', '02'), ('21', '02'), ('01', '20'), ('12', '20')]}
    assert column_arcs == expected_arcs
    print('Local squares: 81 tested, exactly 8 valid; column graph is 2 P3 + 3 K1')

    rows = list(product(range(3), repeat=5))
    adjacency = {s: [t for t in rows if compatible(s, t)] for s in rows}
    arcs = {(s, t) for s in rows for t in adjacency[s]}
    edge_names = [
        ('01010', '20202'), ('20202', '12121'),
        ('10101', '02020'), ('02020', '21212'),
        ('01020', '20212'), ('02010', '21202'),
        ('02021', '10102'), ('02120', '10201'),
        ('02121', '10202'), ('12020', '20101'),
        ('12021', '20102'), ('12120', '20201'),
    ]
    expected_edges = {(digits(a), digits(b)) for a, b in edge_names}
    expected_row_arcs = expected_edges | {(t, s) for s, t in expected_edges}
    assert len(expected_row_arcs) == 24 and arcs == expected_row_arcs
    unseen = set(rows)
    components = []
    while unseen:
        s = min(unseen)
        stack, component = [s], set()
        while stack:
            v = stack.pop()
            if v not in component:
                component.add(v)
                stack.extend(adjacency[v])
        unseen -= component
        components.append(component)
    assert Counter(map(len, components)) == Counter({1: 221, 2: 8, 3: 2})
    assert Counter(map(len, adjacency.values())) == Counter({0: 221, 1: 20, 2: 2})
    # Independent count of four-edge walks in the local column graph.
    col_vertices = list(product(range(3), repeat=2))
    column_counts = {v: 1 for v in col_vertices}
    for _ in range(4):
        column_counts = {t: sum(column_counts[s] for s in col_vertices
                                if (s, t) in column_arcs) for t in col_vertices}
    assert sum(column_counts.values()) == len(arcs) == 24
    print('Ordered row pairs: 59049 tested; 24 arcs, 12 undirected edges')
    print('Exact row graph: 2 P3 + 8 K2 + 221 K1 (243 rows, 22 nonisolated)')
    for a, b in edge_names:
        print('EDGE', a, b)

    examples = [
        ['21202', '02010', '21202', '02010', '21202'],
        ['02010', '21202', '02010', '21202', '02010'],
        ['01010', '20202', '01010', '20202', '12121'],
        ['01010', '20202', '01010', '20202', '01010'],
    ]
    assert all(all(compatible(digits(s), digits(t)) for s, t in zip(x, x[1:]))
               for x in examples)
    print('Fidelity: all four displayed source examples at n=4 satisfy the literal rule')

    data = Path(__file__).with_name('b209723.txt').read_bytes()
    official = {int(n): int(v) for line in data.decode().splitlines()
                if line.strip() and not line.lstrip().startswith('#')
                for n, v in [line.split()]}
    assert sorted(official) == list(range(1, 211))
    max_n = max(official)
    counts = {s: 1 for s in rows}
    raw, a = {}, {0: Fraction(0)}
    for n in range(1, max_n + 1):
        counts = {t: sum(counts[s] for s in adjacency[t]) for t in rows}
        raw[n] = sum(counts.values())
        a[n] = Fraction(raw[n], 4)
        assert a[n].denominator == 1 and a[n] == official[n]
        expected = 3 * 2 ** (n // 2 - 1) + 4 if n % 2 == 0 else 2 ** ((n + 1) // 2) + 4
        assert a[n] == expected
        if n >= 4:
            assert a[n] == a[n - 1] + 2 * a[n - 2] - 2 * a[n - 3]
    direct = [direct_count(n) for n in range(1, 9)]
    assert direct == [raw[n] for n in range(1, 9)]
    print('Independent entry-by-entry counts, n=1..8:', ','.join(map(str, direct)))
    print('a(1)..a(12):', ','.join(str(a[n]) for n in range(1, 13)))
    print('Fidelity: all 210 official b-file pairs match raw array counts / 4')
    print('b209723.txt SHA256:', sha256(data).hexdigest())
    print('Parity formulas: n=1..210; recurrence: n=4..210; all exact')
    # Coefficients of (1-X)(1-2X^2) A = 6X + X^2 - 11X^3.
    denominator = [1, -1, -2, 2]
    numerator = {1: 6, 2: 1, 3: -11}
    for n in range(max_n + 1):
        convolution = sum(Fraction(c) * a[n - j] for j, c in enumerate(denominator) if j <= n)
        assert convolution == numerator.get(n, 0)
    print('Cross-multiplied formal GF: coefficients 0..210 match exactly')
    print('Scope: finite checks only; no all-height proof and no Lean proof probe')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
