#!/usr/bin/env python3
"""Exact checks for Hardin's three determinant-monotone array recurrences.

Run with Python 3; standard library only. Source observations are in source-data.json.
No floating point arithmetic, Lean proofs, or inference from fitted recurrences.
"""
from itertools import product
from pathlib import Path
import json

SPECS = (
    ('A204609', 2, (3, 0, -2), 6, 8, 8),
    ('A204610', 3, (6, -5, -12, 12), 7, 18, 5),
    ('A204611', 4, (9, -8, -72, 96, 128, -128, -64), 10, 34, 4),
)


def determinant_row(u, v):
    return tuple(u[j] * v[j + 1] - u[j + 1] * v[j]
                 for j in range(len(u) - 1))


def nondecreasing(d):
    return all(a <= b for a, b in zip(d, d[1:]))


def transfer(width):
    rows = tuple(product((0, 1), repeat=width))
    initial_pairs = [(u, v, determinant_row(u, v)) for u in rows for v in rows
                     if nondecreasing(determinant_row(u, v))]
    states = sorted({(v, d) for u, v, d in initial_pairs})
    index = {s: i for i, s in enumerate(states)}
    b = [0] * len(states)
    for _, v, d in initial_pairs:
        b[index[v, d]] += 1
    matrix = [[0] * len(states) for _ in states]
    for j, (u, d) in enumerate(states):
        for v in rows:
            e = determinant_row(u, v)
            if nondecreasing(e) and all(a <= b for a, b in zip(d, e)):
                matrix[index[v, e]][j] += 1
    assert all(a in (0, 1) for row in matrix for a in row)
    return states, b, matrix


def matvec(matrix, v):
    return [sum(a * b for a, b in zip(row, v)) for row in matrix]


def brute_counts(width, max_n):
    """Exhaust all integer-coded row words; do not use transfer states or T.

    This second representation uses bit access, rather than tuple rows. Its only
    precomputation is the literal adjacent-minor vector of each ordered row pair.
    """
    r = range(1 << width)
    minors = {}
    for u in r:
        for v in r:
            minors[u, v] = tuple(
                ((u >> (width - j - 1)) & 1) * ((v >> (width - j - 2)) & 1)
                - ((u >> (width - j - 2)) & 1) * ((v >> (width - j - 1)) & 1)
                for j in range(width - 1))
    result = []
    for n in range(1, max_n + 1):
        total = 0
        for word in product(r, repeat=n + 1):
            last = None
            for u, v in zip(word, word[1:]):
                d = minors[u, v]
                if any(d[j] > d[j + 1] for j in range(width - 2)):
                    break
                if last is not None and any(last[j] > d[j] for j in range(width - 1)):
                    break
                last = d
            else:
                total += 1
        result.append(total)
    return result


def polynomial_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def gf_certificate(matrix, b, vectors):
    # q_j are polynomial coefficient lists in ascending degree.
    D = [1, -3, 0, 2]
    expected_P = [16, -1, -16, -8, -4, -1]
    assert polynomial_mul([1, -1], [1, -2, -2]) == D
    q = [[sum(D[t] * vectors[k - t][i]
              for t in range(min(k, len(D) - 1) + 1))
          for k in range(6)] for i in range(len(b))]
    for k in range(7):
        left = [(q[i][k] if k < 6 else 0) for i in range(len(b))]
        if k:
            image = matvec(matrix, [q[i][k - 1] for i in range(len(b))])
            left = [a - c for a, c in zip(left, image)]
        right = [(D[k] if k < len(D) else 0) * v for v in b]
        assert left == right, ('GF vector certificate', k, left, right)
    assert [sum(qi[k] for qi in q) for k in range(6)] == expected_P
    print('GF: D =', D, '; P =', expected_P, '(ascending coefficients)')
    print('GF: q =', q)
    print('GF: (I-xT)q = D*b and sum(q) = P: EXACT_OK')
    # Independently test the cleared scalar identity on every source coefficient.
    for k in range(len(vectors)):
        lhs = sum(D[t] * sum(vectors[k - t])
                  for t in range(min(k, len(D) - 1) + 1))
        assert lhs == (expected_P[k] if k < len(expected_P) else 0)
    print('GF: cleared scalar identity through degree 209: EXACT_OK')


def main():
    sources = json.loads(Path(__file__).with_name('source-data.json').read_text())
    for anum, width, coeffs, threshold, expected_states, brute_max in SPECS:
        source = sources[anum]
        states, b, matrix = transfer(width)
        assert len(states) == expected_states
        assert source['offset'] == 1
        assert source['recurrence_coefficients'] == list(coeffs)
        assert source['strict_threshold'] == threshold
        assert [n for n, _ in source['bfile']] == list(range(1, 211))
        vectors = [b]
        for _ in range(209):
            vectors.append(matvec(matrix, vectors[-1]))
        counts = list(map(sum, vectors))
        assert counts == [a for _, a in source['bfile']], anum
        brute = brute_counts(width, brute_max)
        assert brute == counts[:brute_max], (anum, brute, counts[:brute_max])
        print(anum, 'width =', width, 'state_count =', len(states),
              'edge_count =', sum(map(sum, matrix)))
        print('source revision:', source['revision_line'])
        print('source b-file n=1..210: EXACT_MATCH (210 terms)')
        print('exhaustive direct row words n=1..' + str(brute_max) + ':', brute)
        print('states: index, row, last-minor vector, initial multiplicity')
        for i, ((row, d), bi) in enumerate(zip(states, b)):
            print(i, ''.join(map(str, row)), d, bi)
        print('T columns as lists of target indices (every nonzero entry is 1):')
        for j in range(len(states)):
            print(j, [i for i in range(len(states)) if matrix[i][j]])
        # p(t) = t^d - c1*t^(d-1) - ... - cd.
        # threshold - degree = 3 for all three source recurrences.
        degree = len(coeffs)
        assert threshold - degree == 3
        residual = [vectors[threshold][i] - sum(
            coeffs[k - 1] * vectors[threshold - k][i]
            for k in range(1, degree + 1)) for i in range(len(states))]
        if anum != 'A204611':
            assert residual == [0] * len(states)
            print('certificate: T^3*p(T)*b = 0: EXACT_OK')
        else:
            primitive = [0] * len(states)
            for i in (9, 10, 11):
                primitive[i] = -1
            for i in (29, 32, 33):
                primitive[i] = 1
            assert residual == [20 * a for a in primitive]
            assert matvec(matrix, primitive) == primitive
            assert sum(primitive) == 0
            print('certificate R =', primitive)
            print('certificate: T^3*p(T)*b = 20*R; T*R = R; sum(R) = 0: EXACT_OK')
        print('p coefficients, descending:', [1] + [-a for a in coeffs])
        print('residual vector:', residual)
        for n in range(threshold + 1, 211):
            assert counts[n - 1] == sum(
                coeffs[k - 1] * counts[n - k - 1]
                for k in range(1, degree + 1)), (anum, n)
        print('source recurrence for n=' + str(threshold + 1) + '..210: EXACT_OK')
        if anum == 'A204609':
            gf_certificate(matrix, b, vectors)
        print()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
