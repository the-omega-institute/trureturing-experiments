#!/usr/bin/env python3
"""Exact checks of OEIS A396981, revision 45, and its Conjecture 4.

The table builder uses triangular-block comparisons; the independent recursive
implementation uses the integer-square-root expression in Luschny's source
program. These finite checks do not constitute a proof of the conjecture.
"""

from functools import cache
from math import isqrt


def triangle(r: int) -> int:
    return r * (r + 1) // 2


def make_table(max_n: int, max_q: int) -> tuple[list[list[int]], list[int]]:
    blocks = [0] * (max_q + 1)
    s = 1
    for q in range(1, max_q + 1):
        while triangle(s) < q:
            s += 1
        assert triangle(s - 1) < q <= triangle(s)
        assert 1 <= s <= q
        blocks[q] = s
    rows = [[1] * (max_q + 1)]
    for n in range(1, max_n + 1):
        row = [0] * (max_q + 1)
        for q in range(1, max_q + 1):
            s = blocks[q]
            row[q] = row[q - s] + (n + s - 1) * rows[n - 1][q]
        rows.append(row)
    return rows, blocks


@cache
def source_python(n: int, q: int) -> int:
    # Independent evaluation of the source program, attributed to Peter Luschny.
    if n == 0:
        return 1
    if q <= 0:
        return 0
    s = (isqrt(8 * q) + 1) // 2
    return source_python(n, q - s) + (n + s - 1) * source_python(n - 1, q)


def require_equal(actual: int, expected: int, label: object) -> None:
    if actual != expected:
        raise AssertionError((label, actual, expected))


def main() -> None:
    max_n, max_k, max_m = 14, 50, 5
    max_q = max_k + triangle(max_k + max_m + max_n)
    rows, blocks = make_table(max_n, max_q)

    # Literal published DATA, read by ascending antidiagonals A(d-k,k).
    data = [int(t) for t in (
        "1,0,1,0,1,1,0,2,2,1,0,6,6,3,1,0,24,24,11,4,1,0,120,120,50,18,5,1,0,"
        "720,720,274,96,26,6,1,0,5040,5040,1764,600,154,35,7,1,0,40320,40320,"
        "13068,4320,1044,225,46,8,1,0,362880,362880,109584,35280,8028,1624,326,58,9,1"
    ).split(",")]
    index = 0
    for d in range(10):
        for k in range(d + 1):
            require_equal(rows[d - k][k], data[index], ("DATA", index))
            index += 1
    # The source continues with all 11 entries of antidiagonal 10.
    for k in range(11):
        require_equal(rows[10 - k][k], data[index], ("DATA", index))
        index += 1
    require_equal(index, len(data), "DATA length")
    print(f"FIDELITY published DATA: {index} entries, zero mismatches")

    examples = [
        [1, 1, 1, 1, 1, 1, 1, 1],
        [0, 1, 2, 3, 4, 5, 6, 7],
        [0, 2, 6, 11, 18, 26, 35, 46],
        [0, 6, 24, 50, 96, 154, 225, 326],
        [0, 24, 120, 274, 600, 1044, 1624, 2556],
        [0, 120, 720, 1764, 4320, 8028, 13132, 22212],
        [0, 720, 5040, 13068, 35280, 69264, 118124, 212976],
        [0, 5040, 40320, 109584, 322560, 663696, 1172700, 2239344],
        [0, 40320, 362880, 1026576, 3265920, 6999840, 12753576, 25659360],
    ]
    for n, expected_row in enumerate(examples):
        for q, expected in enumerate(expected_row):
            require_equal(rows[n][q], expected, ("EXAMPLE", n, q))
    print("FIDELITY published rectangular examples: 72 entries, zero mismatches")

    for q in range(1, max_q + 1):
        require_equal(blocks[q], (isqrt(8 * q) + 1) // 2, ("BLOCK", q))
    print(f"BLOCK triangular inequalities versus source isqrt: {max_q} indices, zero mismatches")
    count = 0
    for n in range(9):
        for q in range(81):
            require_equal(rows[n][q], source_python(n, q), ("SOURCE-PYTHON", n, q))
            count += 1
    print(f"FIDELITY independent recursive source evaluation: {count} pairs, zero mismatches")

    count = 0
    for n in range(1, 15):
        for q in range(1, 601):
            require_equal(rows[n][q] - rows[n][q - 1], rows[n - 1][q + blocks[q]],
                          ("DIFFERENCE", n, q))
            count += 1
    print(f"DIFFERENCE LAW n=1..14, q=1..600: {count} pairs, zero mismatches")

    count = 0
    for n in range(max_n + 1):
        for k in range(max_k + 1):
            for m in range(max_m + 1):
                rhs_q = triangle(k + m + 1) - m
                require_equal(rhs_q, triangle(k + m) + k + 1, ("RHS-INDEX", k, m))
                assert rhs_q >= 0
                lhs = sum(rows[n - j][k + triangle(k + m + j)] for j in range(n + 1))
                require_equal(lhs, rows[n][rhs_q], ("CONJECTURE-4", n, k, m))
                count += 1
    print(f"CONJECTURE 4 n=0..14, k=0..50, m=0..5: {count} triples, zero mismatches")
    print(f"TABLE rows=0..{max_n}, columns=0..{max_q}; exact Python integers only")

    # An unrelated unlabelled source formula is inconsistent with NAME.
    # It is neither a premise nor a target of these checks.
    require_equal(rows[1][triangle(2)], 3, "UNUSED-FORMULA lhs")
    require_equal(4 * 3 // 2, 6, "UNUSED-FORMULA c(4,3), one transposition")
    print("SOURCE NOTE unused Stirling FORMULA at n=1,k=3: NAME gives 3; c(4,3)=6")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
