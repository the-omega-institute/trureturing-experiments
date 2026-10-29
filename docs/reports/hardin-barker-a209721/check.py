#!/usr/bin/env python3
"""Exact checks of the source objects in OEIS A209721 (revision 43).

Standard library only. The stored source table is data for fidelity checks;
array counts are computed independently of that table and of the formulas.
"""
from fractions import Fraction
from itertools import product
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def differences(a, b, c, d):
    return (b - a, d - b, c - d, a - c)


def allowed(a, b, c, d):
    return len(set(differences(a, b, c, d))) == 4


def compatible(top, bottom):
    return all(allowed(top[j], top[j + 1], bottom[j], bottom[j + 1])
               for j in range(2))


def direct_count(n):
    # Literal entry-by-entry enumeration, independently of the row transfer.
    total = 0
    for entries in product(range(3), repeat=3 * (n + 1)):
        if all(len({entries[3*i+j+1] - entries[3*i+j],
                    entries[3*(i+1)+j+1] - entries[3*i+j+1],
                    entries[3*(i+1)+j] - entries[3*(i+1)+j+1],
                    entries[3*i+j] - entries[3*(i+1)+j]}) == 4
               for i in range(n) for j in range(2)):
            total += 1
    return total


def main():
    source = {}
    for line in (ROOT / 'b209721.txt').read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            n, term = map(int, line.split())
            source[n] = term
    assert sorted(source) == list(range(1, 211))
    squares = {s for s in product(range(3), repeat=4) if allowed(*s)}
    expected_squares = {(0,1,2,0), (0,2,1,0), (1,0,0,2), (2,0,0,1),
                        (2,1,0,2), (2,0,1,2), (1,2,2,0), (0,2,2,1)}
    assert squares == expected_squares
    assert all(set(differences(*s)) == {-2, -1, 1, 2} for s in squares)
    rows = tuple(product(range(3), repeat=3))
    edges = {(u, v) for u in rows for v in rows if compatible(u, v)}
    expected_undirected = [('010', '202'), ('121', '202'), ('101', '020'),
                           ('020', '212'), ('021', '102'), ('120', '201')]
    row = lambda s: tuple(map(int, s))
    expected_edges = {(row(u), row(v)) for u, v in expected_undirected}
    expected_edges |= {(v, u) for u, v in expected_edges}
    assert edges == expected_edges
    isolated = [u for u in rows if not any((u, v) in edges for v in rows)]
    assert len(isolated) == 17
    counts = [len(rows)]  # walks of length 0, excluded from the source sequence
    current = {u: 1 for u in rows}
    for n in range(1, 211):
        current = {v: sum(current[u] for u in rows if (u, v) in edges) for v in rows}
        counts.append(sum(current.values()))
    a = {n: Fraction(counts[n], 4) for n in source}
    assert all(a[n] == source[n] for n in source)
    assert [a[n] for n in range(1, 9)] == [3, 4, 5, 7, 9, 13, 17, 25]
    direct = {n: direct_count(n) for n in (1, 2, 3)}
    assert all(direct[n] == counts[n] for n in direct)
    for n in range(1, 211):
        closed = (3 * 2**(n//2 - 1) + 1 if n % 2 == 0
                  else 2**((n+1)//2) + 1)
        assert a[n] == closed
        assert a[n].denominator == 1
    for n in range(3, 211):
        assert a[n] == 2 * a[n-2] - 1
    for n in range(4, 211):
        assert a[n] == a[n-1] + 2*a[n-2] - 2*a[n-3]
    # Exact convolution with D = (1-x)(1-2x^2), no rational-function numerics.
    coeff = lambda n: a.get(n, Fraction(0))
    denominator = [1, -1, -2, 2]
    numerator = {1: 3, 2: 1, 3: -5}
    for n in range(211):
        assert sum(Fraction(denominator[j]) * coeff(n-j)
                   for j in range(4)) == numerator.get(n, 0)
    # Check the generic P3 two-step extension property consumed by the argument.
    p3 = {(0,1), (1,0), (1,2), (2,1)}
    assert all(sum((u,v) in p3 and (v,w) in p3
                   for v in range(3) for w in range(3)) == 2 for u in range(3))
    # Polynomial identity used to combine the parity generating series.
    def add(p, q):
        return [sum(v[i] if i < len(v) else 0 for v in (p, q))
                for i in range(max(len(p), len(q)))]
    def mul(p, q):
        return [sum(p[j]*q[i-j] for j in range(len(p)) if 0 <= i-j < len(q))
                for i in range(len(p)+len(q)-1)]
    assert add(mul([0,2,3], [1,-1]), mul([0,1], [1,0,-2])) == [0,3,1,-5]
    print(json.dumps({
        'source': 'OEIS A209721 revision 43, 2022-07-15',
        'source_bfile_sha256': hashlib.sha256((ROOT/'b209721.txt').read_bytes()).hexdigest(),
        'squares_checked': 81,
        'allowed_squares': sorted(map(list, squares)),
        'ordered_row_pairs_checked': 729,
        'oriented_edges': len(edges),
        'undirected_edges': expected_undirected,
        'isolated_rows': [''.join(map(str, u)) for u in isolated],
        'direct_arrays_enumerated': {n: 3**(3*(n+1)) for n in direct},
        'direct_raw_counts': direct,
        'raw_counts_first_eight': counts[1:9],
        'normalized_counts_first_eight': [int(a[n]) for n in range(1,9)],
        'source_bfile_terms_matched': len(source),
        'recurrence_checked': [4,210],
        'gf_coefficients_checked': [0,210],
        'parity_formula_checked': [1,210],
        'p3_two_step_extensions_per_start': [2,2,2],
        'scope': 'Exact finite checks; no all-height or Lean proof is claimed.'
    }, indent=2, sort_keys=True))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
