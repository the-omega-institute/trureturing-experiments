#!/usr/bin/env python3
"""Exact finite checks for Kurkov's A325804 gap conjecture; not a proof."""

import argparse
from collections import Counter
from hashlib import sha256
from pathlib import Path


# OEIS A325804 revision 45, fields S/T/U, offset 1.
SOURCE_SUPPORT = [
    0, 1, 3, 6, 7, 12, 14, 15, 24, 25, 28, 29, 30, 31, 48, 50, 51, 56,
    57, 58, 60, 61, 62, 63, 96, 97, 100, 101, 102, 103, 112, 113, 114,
    115, 116, 117, 120, 121, 122, 123, 124, 125, 126, 127, 192, 194,
    195, 200, 201, 202, 204, 205, 206, 207, 224, 225, 226, 228,
]
# Table in Drexler-Lemire--Shallit, arXiv:1402.3091v1, Introduction.
SOURCE_NOTES = [0, 1, -1, 2, 1, 0, -2, 3, -1, 2, 0, 1, 2, -1, -3, 4,
                1, 0, -2, 3]
# OEIS A094373, fields S/T/U, offset 0.
SOURCE_GAPS = [
    1, 2, 3, 5, 9, 17, 33, 65, 129, 257, 513, 1025, 2049, 4097,
    8193, 16385, 32769, 65537, 131073, 262145, 524289, 1048577,
    2097153, 4194305, 8388609, 16777217, 33554433, 67108865,
    134217729, 268435457, 536870913, 1073741825, 2147483649,
    4294967297, 8589934593,
]


def require(test, message):
    if not test:
        raise AssertionError(message)


def bit_fold(n):
    """Independent MSB-first evaluation and all-prefix legality."""
    q = 0
    legal = True
    for bit in bin(n)[2:] if n else "":
        q = q + 1 if bit == "1" else -q
        legal = legal and q != -1
    return q, legal


def literal_product(n, notes):
    """Literal source product, with P(0)=1; bit_length avoids real logs."""
    value = 1
    for k in range(n.bit_length()):
        value *= 1 + notes[n // (2 ** k)]
    return value


def legal_suffixes(q, length):
    """Exhaust all words; this does not use the proposed extremum formula."""
    states = [(0, q)]
    for _ in range(length):
        new = []
        for value, state in states:
            if -state != -1:
                new.append((2 * value, -state))
            if state + 1 != -1:
                new.append((2 * value + 1, state + 1))
        states = new
    return [value for value, _ in states]


def proposed_extrema(q, m):
    if m == 0:
        return 0, 0
    if q == 0:
        return 2 ** m - 1, 2 ** (m - 1)
    if q <= -3 or m <= q - 2:
        return 2 ** m - 1, 0
    require(q >= 2 and m >= q - 1, "missing state case")
    return 2 ** m - 1 - 2 ** (m - q + 1), 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bfile", type=Path,
                        help="optional downloaded https://oeis.org/A325804/b325804.txt")
    args = parser.parse_args()
    bound = 2 ** 20
    notes = [0] * bound
    products = [1] * bound
    for n in range(1, bound):
        notes[n] = notes[n // 2] + 1 if n % 2 else -notes[n // 2]
        products[n] = (1 + notes[n]) * products[n // 2]
    support = [n for n in range(bound) if products[n] != 0]
    require(notes[:len(SOURCE_NOTES)] == SOURCE_NOTES, "source note table")
    require(support[:len(SOURCE_SUPPORT)] == SOURCE_SUPPORT, "source support table")
    gaps = [1] + [2 ** e + 1 for e in range(len(SOURCE_GAPS) - 1)]
    require(gaps == SOURCE_GAPS, "source A094373 table")
    print("FIDELITY source_notes=20 source_support_terms=58 source_gap_terms=35")
    for n in range(4096):
        require(literal_product(n, notes) == products[n], f"literal P at {n}")
        q, legal = bit_fold(n)
        require(q == notes[n], f"independent note fold at {n}")
        require(legal == (products[n] != 0), f"independent prefix legality at {n}")
    print("INDEPENDENT_FOLD_AND_LITERAL_PRODUCT n=0..4095 cases=4096")
    histogram = Counter()
    for x, y in zip(support, support[1:]):
        d = y - x
        require(d > 0, "increasing enumeration")
        require(d == 1 or (d > 1 and (d - 1) & (d - 2) == 0),
                f"conjecture failure: {x}, {y}, gap={d}")
        histogram[d] += 1
    require(len(support) == 191674, "candidate support count")
    print(f"GAPS n=0..{bound - 1} support={len(support)} adjacent_pairs={len(support) - 1}")
    print("GAP_HISTOGRAM " + " ".join(f"{d}:{histogram[d]}" for d in sorted(histogram)))
    cases = 0
    for q in range(-9, 10):
        if q in (-2, -1, 1):
            continue
        for m in range(10):
            left = legal_suffixes(-q, m)
            right = legal_suffixes(q + 1, m)
            require(bool(left) and bool(right), f"empty branch q={q}, m={m}")
            require((max(left), min(right)) == proposed_extrema(q, m),
                    f"extremum formula q={q}, m={m}")
            cases += 1
    require(cases == 160, "suffix case count")
    print("SUFFIX_EXTREMA q=-9..9 excluding -2,-1,1 lengths=0..9 cases=160")
    if args.bfile is not None:
        data = args.bfile.read_bytes()
        rows = []
        for line in data.decode().splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                i, value = map(int, line.split())
                rows.append((i, value))
        require([i for i, _ in rows] == list(range(1, len(rows) + 1)),
                "b-file index coverage")
        require([value for _, value in rows] == support[:len(rows)], "b-file fidelity")
        require(rows[-1] == (13496, 65535), "candidate b-file last row")
        print(f"BFILE matched_rows={len(rows)} last_index={rows[-1][0]} last_value={rows[-1][1]}")
        print("BFILE_SHA256 " + sha256(data).hexdigest())
    print("KERNEL_CHECKED false; finite checks do not establish the universal claim")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
