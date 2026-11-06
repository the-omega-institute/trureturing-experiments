#!/usr/bin/env python3
"""Exact preregistration checks for OEIS A185761--A185764; no Lean proof."""
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
# a(n) = sum(coeff[j-1] * a(n-j), j=1..d), starting at start.
# The final two lists are ascending coefficients of the published OGF numerator/denominator.
CASES = {
    'A185761': (2, [1, 4], 3, [0, 16, 24], [1, -1, -4]),
    'A185762': (3, [1, 6, 1], 4, [0, 40, 80, 8], [1, -1, -6, -1]),
    'A185763': (4, [1, 10, 3, -5, 1], 7, [0, 104, 264, 0, -128, 56, -8],
                [1, -1, -10, -3, 5, -1]),
    'A185764': (5, [1, 17, 9, -22, -1, 1, -6], 10,
                [0, 264, 864, -344, -1992, -552, 40, 344, 648, 144],
                [1, -1, -17, -9, 22, 1, -1, 6]),
}


def bfile(name):
    rows = [tuple(map(int, line.split())) for line in (ROOT / name).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith('#')]
    assert all(len(row) == 2 for row in rows)
    assert [n for n, _ in rows] == list(range(1, len(rows) + 1))
    return dict(rows)


def distinct_neighbors(t):
    return all(x != y for x, y in zip(t, t[1:]))


def traces(a, b):
    return tuple(a[i] + b[i + 1] for i in range(len(a) - 1))


def complement(s):
    b, t = s
    return tuple(1 - x for x in b), tuple(2 - x for x in t)


def canonical(s):
    return min(s, complement(s))


def next_states(s, columns):
    b, t = s
    for c in columns:
        u = tuple(b[i] + c[i + 1] for i in range(len(b)))
        if distinct_neighbors(u) and all(x != y for x, y in zip(t, u)):
            yield c[:-1], u


def transfer(width):
    columns = list(product(range(2), repeat=width))
    initial = Counter()
    for a, b in product(columns, repeat=2):
        t = traces(a, b)
        if distinct_neighbors(t):
            initial[(b[:-1], t)] += 1
    raw = sorted(initial)
    states = sorted({canonical(s) for s in raw})
    index = {s: i for i, s in enumerate(states)}
    matrix = []
    for s in states:
        transitions = Counter(canonical(t) for t in next_states(s, columns))
        matrix.append([transitions[t] for t in states])
    # Check closure and every representative, rather than assuming lumpability.
    for s in raw:
        assert complement(complement(s)) == s and complement(s) != s
        assert complement(s) in initial and initial[s] == initial[complement(s)]
        successors = list(next_states(s, columns))
        assert all(t in initial for t in successors)
        row = Counter(canonical(t) for t in successors)
        assert [row[t] for t in states] == matrix[index[canonical(s)]]
    vector = [0] * len(states)
    for s, multiplicity in initial.items():
        vector[index[canonical(s)]] += multiplicity
    return len(raw), states, matrix, vector


def matvec(matrix, vector):
    return [sum(x * y for x, y in zip(row, vector)) for row in matrix]


def rowmul(vector, matrix):
    return [sum(vector[i] * matrix[i][j] for i in range(len(vector)))
            for j in range(len(vector))]


def transfer_counts(matrix, initial, limit):
    v = list(initial)
    values = {}
    for n in range(1, limit + 1):
        values[n] = sum(v)
        v = rowmul(v, matrix)
    return values


def full_pair_counts(width, limit):
    """Separate uncompressed transfer retaining both entire preceding columns."""
    columns = list(product(range(2), repeat=width))
    states = [(a, b) for a, b in product(columns, repeat=2)
              if distinct_neighbors(traces(a, b))]
    index = {s: i for i, s in enumerate(states)}
    successors = []
    for a, b in states:
        old = traces(a, b)
        successors.append([index[(b, c)] for c in columns
                           if (b, c) in index
                           and all(x != y for x, y in zip(old, traces(b, c)))])
    v = [1] * len(states)
    values = {}
    for n in range(1, limit + 1):
        values[n] = sum(v)
        w = [0] * len(states)
        for i, js in enumerate(successors):
            for j in js:
                w[j] += v[i]
        v = w
    return len(states), values


def literal_count(n, width):
    """Enumerate original (n+1)-by-width arrays, with no transfer representation."""
    good = 0
    for bits in product(range(2), repeat=(n + 1) * width):
        t = [[bits[i * width + j] + bits[(i + 1) * width + j + 1]
              for j in range(width - 1)] for i in range(n)]
        if (all(t[i][j] != t[i + 1][j] for i in range(n - 1) for j in range(width - 1))
                and all(t[i][j] != t[i][j + 1] for i in range(n) for j in range(width - 2))):
            good += 1
    return good


def polynomial_certificate(matrix, coefficients):
    """Find b and w=T^b*1 with q(T)w=0, q monic with the given recurrence."""
    d = len(coefficients)
    q = [-coefficients[d - 1 - j] for j in range(d)] + [1]
    w = [1] * len(matrix)
    for burn in range(len(matrix) + 1):
        powers = [w]
        for _ in range(d):
            powers.append(matvec(matrix, powers[-1]))
        residual = [sum(q[j] * powers[j][i] for j in range(d + 1))
                    for i in range(len(matrix))]
        if not any(residual):
            return q, burn, w, residual
        w = matvec(matrix, w)
    raise AssertionError('No exact vector certificate')


def quadratic_mul(a, b):
    return a[0] * b[0] + 17 * a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def barker_closed(n):
    plus = (Fraction(1), Fraction(0))
    minus = plus
    for _ in range(n):
        plus = quadratic_mul(plus, (Fraction(1, 2), Fraction(1, 2)))
        minus = quadratic_mul(minus, (Fraction(1, 2), Fraction(-1, 2)))
    a = quadratic_mul((Fraction(13), Fraction(3)), plus)
    b = quadratic_mul((Fraction(-13), Fraction(3)), minus)
    # Divide (a+b) by sqrt(17), in Q[sqrt(17)].
    return a[1] + b[1], (a[0] + b[0]) / 17


def main():
    entries = json.loads((ROOT / 'oeis-entries.json').read_text())
    all_counts = {}
    for anum, (width, coeffs, start, numerator, denominator) in CASES.items():
        published = bfile('b' + anum[1:] + '.txt')
        assert len(published) == 200
        raw, states, matrix, initial = transfer(width)
        assert len(states) <= 64
        assert (raw, len(states)) == {2: (6, 3), 3: (16, 8), 4: (40, 20), 5: (104, 52)}[width]
        limit = max(published)
        counted = transfer_counts(matrix, initial, limit)
        pair_size, independent = full_pair_counts(width, limit)
        assert counted == independent == published
        source_prefix = list(map(int, entries[anum]['data'].split(',')))
        assert [counted[i] for i in range(1, len(source_prefix) + 1)] == source_prefix
        brute = {}
        for n in range(1, 8):
            if (n + 1) * width <= 16:
                brute[n] = literal_count(n, width)
                assert brute[n] == counted[n]
        for n in range(start, limit + 1):
            assert counted[n] == sum(c * counted[n - j] for j, c in enumerate(coeffs, 1))
        series = {0: 0, **counted}
        for n in range(limit + 1):
            value = sum(c * series[n - j] for j, c in enumerate(denominator) if j <= n)
            assert value == (numerator[n] if n < len(numerator) else 0)
        q, burn, w, residual = polynomial_certificate(matrix, coeffs)
        eventual_start = len(coeffs) + burn + 1
        for n in range(start, eventual_start):
            assert counted[n] == sum(c * counted[n - j] for j, c in enumerate(coeffs, 1))
        if anum == 'A185761':
            for n in published:
                assert barker_closed(n) == (Fraction(counted[n]), Fraction(0))
        matrix_bytes = json.dumps(matrix, separators=(',', ':')).encode()
        print(f'{anum}: width={width}; raw_compressed_states={raw}; complement_orbit_states={len(states)}; full_pair_states={pair_size}')
        print(f'  source_revision={entries[anum]["revision"]}; source_updated={entries[anum]["time"]}')
        print(f'  all_bfile_terms=200; range=1..200; source_prefix_terms={len(source_prefix)}; fidelity_and_independent_transfer=OK')
        print(f'  literal_bruteforce={brute}')
        print(f'  recurrence_coefficients={coeffs}; first_claimed_n={start}; finite_bridge_n={list(range(start, eventual_start))}')
        print(f'  q_coefficients_ascending={q}; burn={burn}; q(T)T^burn*1={residual}')
        print(f'  certificate_w=T^burn*1={w}')
        print(f'  initial_orbit_multiplicities={initial}')
        print(f'  canonical_matrix_sha256={hashlib.sha256(matrix_bytes).hexdigest()}')
        print(f'  OGF_numerator={numerator}; OGF_denominator={denominator}; coefficient_checks=0..200 OK')
        if anum == 'A185761':
            print('  Barker_closed_form_in_Q_sqrt17=OK for every n=1..200')
        all_counts[width] = counted
    # Parent table uses antidiagonals (n,k)=(d-1,1),...,(1,d-1).
    table = bfile('b185769.txt')
    comparisons = 0
    idx, diagonal = 1, 2
    while idx <= max(table):
        for k in range(1, diagonal):
            if idx > max(table):
                break
            n = diagonal - k
            if k + 1 in all_counts:
                assert table[idx] == all_counts[k + 1][n]
                comparisons += 1
            idx += 1
        diagonal += 1
    print(f'A185769 parent table: {comparisons} applicable cells checked among all {len(table)} bfile entries')
    raw6, states6, _, _ = transfer(6)
    assert (raw6, len(states6)) == (264, 132)
    print('Width 6 excluded: raw_compressed_states=264; complement_orbit_states=132 > 64; first four smallest widths selected')
    print('Evidence is exact integer/rational arithmetic; no floating point and no Lean proof.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
