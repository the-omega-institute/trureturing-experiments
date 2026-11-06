#!/usr/bin/env python3
"""Exact checks for Mathar's whole generating-function conjecture, OEIS A282523.

Python 3.10+ standard library only. No floating point, external prover, or network.
The literal-cell predicate and the row predicate are implemented independently.
This is computational evidence and finite certificate data, not a Lean proof.
Source: https://oeis.org/A282523/internal, revision 9; OEIS data CC BY-SA 4.0.
The proposed certificate was supplied by research search seat codex-search-b.
"""

from collections import Counter


CLASS = [
    [0, 1, 2, 3, 1, 4, 3, 5],
    [0, 6, 7, 5, 1, 8, 5, 9],
    [0, 6, 7, 5, 6, 10, 5, 9],
    [0, 11, 5, 9, 6, 12, 9, 9],
    [0, 1, 7, 5, 6, 8, 5, 9],
    [0, 6, 5, 9, 6, 10, 9, 9],
    [0, 6, 5, 9, 11, 12, 9, 9],
    [0, 11, 9, 9, 11, 5, 9, 9],
]
M = [
    [1, 2, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 0, 2, 1, 1, 1, 1, 0, 0, 0],
    [1, 0, 0, 0, 0, 2, 2, 1, 0, 0, 1, 0, 0],
    [1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0, 1, 2, 0, 0, 3, 1, 0, 0],
    [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 1, 2, 0, 0, 1, 1, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 1, 0, 0],
    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
]
F = [1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1]
ALPHA = [1, 2, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0]
Q = [1, -5, -2, 1, 24, -13, 0, 2, -6]
N = [
    [1, 3, 1, -24, 13, 0, -2, 6],
    [1, 2, -5, -1, 3, 0, 0, 0],
    [1, 2, -6, 3, 2, -2, 0, 0],
    [1, -1, -5, 5, -3, 0, 3, 0],
    [1, 0, -1, 0, -2, 0, 0, 0],
    [1, -4, 1, 2, 0, 0, 0, 0],
    [1, 1, -1, -2, 1, 0, 0, 0],
    [1, -1, 2, 2, -2, 0, 0, 0],
    [1, 0, -1, 0, -2, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 0, 0],
    [1, 0, -1, 0, -2, 0, 0, 0],
    [1, -3, 3, -3, -1, 3, 0, 0],
    [1, -3, 2, 1, -2, 1, 0, 0],
]
P = [0, 8, 3, -25, -11, 13, -2, 4, 6]
SOURCE_TERMS = [
    8, 43, 206, 1097, 5675, 29433, 153037, 794716, 4128244,
    21444844, 111394775, 578646848, 3005807582, 15613798271,
    81106596833, 421311892656, 2188523765507, 11368386235590,
    59053598944793, 306756604030647, 1593461123288541,
]


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def bit(row, j):
    return (row >> j) & 1 if 0 <= j < 3 else 0


def row_allowed(u, v, w):
    """Check only the occupied cells of the middle row."""
    for j in range(3):
        if bit(v, j):
            degree = bit(v, j - 1) + bit(v, j + 1)
            degree += sum(bit(r, k) for r in (u, w)
                          for k in range(j - 1, j + 2))
            if degree > 2:
                return False
    return True


def literal_neighborhoods(height):
    """Direct Chebyshev-distance definition, with no padded rows or row triples."""
    cells = [(i, j) for i in range(height) for j in range(3)]
    return [sum(1 << q for q, (k, l) in enumerate(cells)
                if (i, j) != (k, l) and max(abs(i - k), abs(j - l)) == 1)
            for i, j in cells]


def literal_valid(board, neighbors, selected=None):
    cells = range(len(neighbors)) if selected is None else selected
    return all(not (board & (1 << i)) or (board & neighbors[i]).bit_count() <= 2
               for i in cells)


def literal_count(height):
    neighbors = literal_neighborhoods(height)
    return sum(literal_valid(board, neighbors) for board in range(1 << (3 * height)))


def matvec(matrix, vector):
    return [sum(a * b for a, b in zip(row, vector)) for row in matrix]


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def main():
    require(len(CLASS) == 8 and all(len(r) == 8 for r in CLASS), 'class shape')
    require(set(sum(CLASS, [])) == set(range(13)), 'all classes represented')
    require(len(M) == len(N) == len(F) == len(ALPHA) == 13, 'certificate shape')
    require(all(len(r) == 13 for r in M), 'matrix shape')
    require(all(len(r) == 8 for r in N), 'polynomial shape')
    require(len(P) == len(Q) == 9 and Q[0] == 1, 'unit denominator')

    triple_neighbors = literal_neighborhoods(3)
    for u in range(8):
        for v in range(8):
            s = CLASS[u][v]
            require(int(row_allowed(u, v, 0)) == F[s], ('terminal', u, v))
            counts = Counter()
            for w in range(8):
                board = u | (v << 3) | (w << 6)
                allowed = row_allowed(u, v, w)
                require(allowed == literal_valid(board, triple_neighbors, range(3, 6)),
                        ('independent local triple', u, v, w))
                if allowed:
                    counts[CLASS[v][w]] += 1
            require([counts[t] for t in range(13)] == M[s], ('quotient', u, v))
    initial = Counter(CLASS[0][v] for v in range(8))
    require([initial[t] for t in range(13)] == ALPHA, 'initial multiplicities')
    print('Local fidelity: 512 triples, 64 terminal checks, 832 quotient entries, 13 initial entries OK')

    for s in range(13):
        for degree in range(9):
            lhs = N[s][degree] if degree < 8 else 0
            if degree:
                lhs -= sum(M[s][t] * N[t][degree - 1] for t in range(13))
            require(lhs == Q[degree] * F[s], ('(I-xM)N=Qf', s, degree))
    for degree in range(9):
        lhs = 0 if degree == 0 else sum(ALPHA[s] * N[s][degree - 1]
                                             for s in range(13))
        require(lhs == P[degree], ('x alpha N=P', degree))
    print('Exact polynomial certificate: 117 column coefficients and 9 numerator coefficients OK')

    raw = [[0] * 64 for _ in range(64)]
    terminal = []
    for u in range(8):
        for v in range(8):
            terminal.append(int(row_allowed(u, v, 0)))
            for w in range(8):
                raw[8 * u + v][8 * v + w] = int(row_allowed(u, v, w))
    initial_raw = [1] * 8 + [0] * 56
    powers = [terminal]
    for _ in range(8):
        powers.append(matvec(raw, powers[-1]))
    require(all(sum(Q[8 - k] * powers[k][s] for k in range(9)) == 0
                for s in range(64)), 'raw degree-eight annihilator')
    print('Independent raw transfer: p(T)f=0 in all 64 coordinates OK')

    raw_vector, quotient_vector = terminal, F
    counts = [0]
    for height in range(1, 211):
        a_raw = dot(initial_raw, raw_vector)
        a_quotient = dot(ALPHA, quotient_vector)
        require(a_raw == a_quotient, ('raw/quotient', height))
        counts.append(a_raw)
        raw_vector = matvec(raw, raw_vector)
        quotient_vector = matvec(M, quotient_vector)
    require(counts[1:22] == SOURCE_TERMS, 'OEIS revision 9 displayed terms')
    print('Source fidelity a(1..21):', ','.join(map(str, counts[1:22])))
    for height in range(1, 7):
        literal = literal_count(height)
        require(literal == counts[height], ('literal-cell enumeration', height))
        print(f'Literal exhaustive count: n={height}, boards={1 << (3 * height)}, a(n)={literal}')
    for degree in range(211):
        coefficient = sum(Q[k] * counts[degree - k] for k in range(min(8, degree) + 1))
        require(coefficient == (P[degree] if degree < 9 else 0),
                ('Q times generating series coefficient', degree))
    print('Raw/quotient equality through n=210 and Q*A=P through degree 210 OK')
    print('Boundary: finite exact evidence only; no kernel verification of the all-height counting bridge')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
