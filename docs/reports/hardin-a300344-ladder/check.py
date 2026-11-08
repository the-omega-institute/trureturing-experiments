#!/usr/bin/env python3
"""Exact checks for Hardin--Plouffe OEIS A300344; Python 3 standard library.

Source: https://oeis.org/A300344, revision 12.  This is a finite computational
cross-check, not a Lean proof.  Arrays are labelled, with no periodic boundary.
"""
from itertools import product

PUBLISHED = [1, 2, 6, 15, 37, 90, 223, 550, 1355, 3341, 8237, 20307,
             50064, 123426, 304290, 750185, 1849478, 4559634, 11241151,
             27713513, 68323858, 168443083, 415273274, 1023799191,
             2524036217, 6222664445, 15341124083, 37821433280, 93243546406]
EXAMPLE_ROWS = [
    "00 01 00 01 00 00 01 00 00 01",
    "01 01 10 01 11 01 01 10 11 01",
    "10 10 10 01 01 10 01 11 01 01",
    "01 10 10 01 01 01 10 01 00 00",
    "11 11 11 00 01 01 10 01 11 11",
]
EXPECTED_N = [
    [0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0, 1, 1],
    [0, 0, 0, 0, 1, 1, 0, 0],
    [0, 1, 1, 1, 0, 0, 0, 0],
    [0, 0, 1, 1, 0, 0, 0, 0],
    [0, 0, 0, 0, 1, 1, 1, 0],
    [0, 0, 0, 0, 0, 1, 1, 1],
    [0, 0, 0, 1, 0, 0, 0, 0],
]
EXPECTED_U = [0, 1, 1, 1, 0, 1, 0, 0]
EXPECTED_V = [0, 0, 0, 1, 1, 1, 0, 1]
EXPECTED_CONTINUATIONS = [
    [0, 0, 0, 1, 1, 1, 0, 1],
    [0, 1, 2, 1, 1, 2, 2, 1],
    [0, 3, 3, 4, 3, 5, 5, 1],
    [0, 6, 8, 10, 7, 13, 11, 4],
    [0, 15, 20, 24, 18, 31, 28, 10],
    [0, 38, 49, 59, 44, 77, 69, 24],
]


def neighbor_lists(n, antidiagonal=True):
    """Literal coordinate test, independent of the row transfer predicate."""
    cells = list(product(range(n), range(2)))
    neighbors = []
    for i, j in cells:
        adjacent = []
        for r, c in cells:
            horizontal = i == r and abs(j - c) == 1
            vertical = j == c and abs(i - r) == 1
            diagonal = ((i + 1 == r and c + 1 == j) or
                        (r + 1 == i and j + 1 == c)) if antidiagonal else (
                            (i + 1 == r and j + 1 == c) or
                            (r + 1 == i and c + 1 == j))
            if horizontal or vertical or diagonal:
                adjacent.append(2 * r + c)
        assert 2 * i + j not in adjacent
        neighbors.append(tuple(adjacent))
    assert all(p in neighbors[q] for p, qs in enumerate(neighbors) for q in qs)
    return neighbors


def literal_valid(mask, neighbors):
    """Bit k is the entry at the k-th row-major cell; first bit must be zero."""
    if mask & 1:
        return False
    for p, qs in enumerate(neighbors):
        equal = sum(((mask >> p) & 1) == ((mask >> q) & 1) for q in qs)
        if equal not in (1, 2):
            return False
    return True


def bits(r):
    return r // 2, r % 2  # 0=00, 1=01, 2=10, 3=11


def local(a, b, c):
    """Two equal-neighbor degrees in the middle row, omitting absent rows."""
    b0, b1 = bits(b)
    d0 = int(b0 == b1)
    d1 = int(b0 == b1)
    if a is not None:
        a0, a1 = bits(a)
        d0 += int(b0 == a0) + int(b0 == a1)
        d1 += int(b1 == a1)
    if c is not None:
        c0, c1 = bits(c)
        d0 += int(b0 == c0)
        d1 += int(b1 == c0) + int(b1 == c1)
    return d0 in (1, 2) and d1 in (1, 2)


def row_valid(mask, n):
    rows = [2 * ((mask >> (2 * i)) & 1) + ((mask >> (2 * i + 1)) & 1)
            for i in range(n)]
    return not (mask & 1) and all(
        local(rows[i - 1] if i else None, rows[i],
              rows[i + 1] if i + 1 < n else None) for i in range(n))


def reflected_normalized(mask, n):
    """Reflect columns, then complement if necessary to fix upper left zero."""
    reflected = sum(((mask >> (2 * i + 1 - j)) & 1) << (2 * i + j)
                    for i in range(n) for j in range(2))
    return reflected ^ ((1 << (2 * n)) - 1) if reflected & 1 else reflected


def matvec(m, v):
    return [sum(x * y for x, y in zip(row, v)) for row in m]


def rowvec(v, m):
    return [sum(v[i] * m[i][j] for i in range(len(v))) for j in range(len(v))]


def model(quotient):
    states = list(product(range(2 if quotient else 4), range(4)))
    index = {s: i for i, s in enumerate(states)}
    m = [[0] * len(states) for _ in states]
    for i, (a, b) in enumerate(states):
        for c in range(4):
            if local(a, b, c):
                nxt = (3 - b, 3 - c) if quotient and b >= 2 else (b, c)
                m[i][index[nxt]] += 1  # Count actual extensions, including multiplicity.
    u = [int(a < 2 and local(None, a, b)) for a, b in states]
    v = [int(local(a, b, None)) for a, b in states]
    return m, u, v


def transfer_counts(m, u, v, through):
    counts = [0, 1]  # Empty array excluded; only row 00 at n=1.
    row = u[:]
    for n in range(2, through + 1):
        counts.append(sum(x * y for x, y in zip(row, v)))
        row = rowvec(row, m)
    return counts


def main():
    # Check the row rule against literal degrees for each possible local triple.
    ns3 = neighbor_lists(3)
    for a, b, c in product(range(4), repeat=3):
        values = bits(a) + bits(b) + bits(c)
        middle = [sum(values[p] == values[q] for q in ns3[p]) for p in (2, 3)]
        assert local(a, b, c) == all(d in (1, 2) for d in middle)
        assert local(a, b, c) == local(3 - a, 3 - b, 3 - c)
    print("LOCAL_RULE: 64 triples agree with literal coordinate neighbors and complement")

    native_counts, main_counts = [], []
    for n in range(1, 10):
        anti, main_diag = neighbor_lists(n), neighbor_lists(n, False)
        valid_anti, valid_main = set(), set()
        for mask in range(0, 1 << (2 * n), 2):
            ok = literal_valid(mask, anti)
            assert ok == row_valid(mask, n), (n, mask)
            if ok:
                valid_anti.add(mask)
            if literal_valid(mask, main_diag):
                valid_main.add(mask)
        assert {reflected_normalized(mask, n) for mask in valid_anti} == valid_main
        assert all(reflected_normalized(reflected_normalized(mask, n), n) == mask
                   for mask in valid_anti)
        native_counts.append(len(valid_anti))
        main_counts.append(len(valid_main))
        print(f"LITERAL n={n}: antidiagonal={len(valid_anti)}, main_diagonal={len(valid_main)}")
    assert native_counts == PUBLISHED[:9]
    assert main_counts == native_counts
    print("ORIENTATION: the two counts coincide; source examples distinguish the rules")

    examples = [sum(int(bit) << (2 * i + j)
                    for i, row in enumerate(EXAMPLE_ROWS)
                    for j, bit in enumerate(row.split()[k])) for k in range(10)]
    anti_examples = [literal_valid(mask, neighbor_lists(5)) for mask in examples]
    main_examples = [literal_valid(mask, neighbor_lists(5, False)) for mask in examples]
    assert all(anti_examples)
    assert main_examples == [False, False, False, False, False, False, True, False, True, False]
    print("SOURCE_EXAMPLES: antidiagonal=10/10, main_diagonal=2/10")

    n, u, v = model(True)
    assert (n, u, v) == (EXPECTED_N, EXPECTED_U, EXPECTED_V)
    print("QUOTIENT_STATES:", list(product(range(2), range(4))))
    print("N:", n)
    print("u:", u)
    print("v:", v)
    columns, col = [], v[:]
    for _ in range(6):
        columns.append(col)
        col = matvec(n, col)
    assert columns == EXPECTED_CONTINUATIONS
    for k, col in enumerate(columns):
        print(f"N^{k} v:", col)
    residual = [columns[5][i] - columns[4][i] - 2 * columns[3][i] -
                3 * columns[2][i] - 2 * columns[1][i] - columns[0][i]
                for i in range(8)]
    assert residual == [0] * 8
    print("ORDER_FIVE_RESIDUAL:", residual)

    qcounts = transfer_counts(n, u, v, 210)
    m16, u16, v16 = model(False)
    counts16 = transfer_counts(m16, u16, v16, 210)
    assert qcounts == counts16
    assert qcounts[1:10] == native_counts
    assert qcounts[1:30] == PUBLISHED
    for k in range(6, 211):
        assert qcounts[k] == (qcounts[k - 1] + 2 * qcounts[k - 2] +
                              3 * qcounts[k - 3] + 2 * qcounts[k - 4] + qcounts[k - 5])
    denominator = [1, -1, -2, -3, -2, -1]
    numerator = [0, 1, 1, 2, 2, 2]
    for k in range(211):
        coeff = sum(c * qcounts[k - i] for i, c in enumerate(denominator) if i <= k)
        assert coeff == (numerator[k] if k < len(numerator) else 0)
    print("SOURCE_TERMS: all 29 displayed terms agree with the independent 16- and 8-state models")
    print("TRANSFER_CHECK: 16 and 8 states agree for n=1..210")
    print("RECURRENCE_CHECK: n=6..210; initial values", qcounts[1:7])
    print("GENERATING_FUNCTION_CHECK: coefficients 0..210")
    print("EVIDENCE_SCOPE: finite exact arithmetic and certificate checks; no Lean proof")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
