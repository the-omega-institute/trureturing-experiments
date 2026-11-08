#!/usr/bin/env python3
"""Exact independent checks for the A225620 COMMENTS recurrence.

The source fixture is the complete DATA field of OEIS A225620, revision 98.
No network access, third-party packages, or floating-point arithmetic is used.
"""

from functools import cache
from itertools import accumulate

OEIS_REVISION = 98
OEIS_TERMS = (1, 2, 3, 4, 6, 7, 8, 10, 12, 14, 15, 16, 20, 24, 26, 28, 30, 31, 32, 36, 40, 42, 48, 52, 56, 58, 60, 62, 63, 64, 72, 80, 84, 96, 100, 104, 106, 112, 116, 120, 122, 124, 126, 127, 128, 136, 144, 160, 164, 168, 170, 192, 200, 208, 212, 224, 228, 232, 234, 240, 244, 248, 250, 252, 254, 255)
LIMIT = 1 << 15
CODE_BOUND = 1 << 12
PARTITION_COUNTS = (1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176)


def partitions(total, ceiling=None):
    """Independently enumerate decreasing positive partitions of total."""
    if total == 0:
        yield ()
        return
    if ceiling is None:
        ceiling = total
    for first in range(min(total, ceiling), 0, -1):
        for rest in partitions(total - first, first):
            yield (first,) + rest


def source_code(s):
    """Literal source instructions: blocks 0^(j-1)1, concatenate, mirror."""
    assert s and all(x > 0 for x in s)
    blocks = ["0" * (j - 1) + "1" for j in s]
    return int("".join(blocks)[::-1], 2)


def cumulative_code(s):
    return sum(1 << (partial_sum - 1) for partial_sum in accumulate(s))


def decode_code(code):
    positions = [i for i in range(code.bit_length()) if code & (1 << i)]
    assert positions
    return (positions[0] + 1,) + tuple(y - x for x, y in zip(positions, positions[1:]))


def staircase(n):
    exponents = [i for i in range(n.bit_length() - 1, -1, -1) if n & (1 << i)]
    count = len(exponents)
    return tuple(e - count + i + 2 for i, e in enumerate(exponents))


def inverse_staircase(s):
    count = len(s)
    return sum(1 << (part + count - i - 2) for i, part in enumerate(s))


@cache
def b(n):
    assert n > 0
    highest = 1 << (n.bit_length() - 1)
    if n == highest:
        return n
    remainder = n - highest
    zero_count = bin(n)[2:].count("0")
    assert remainder > 0
    return (1 << zero_count) * (2 * b(remainder) + 1)


def main():
    assert source_code((2, 2, 1, 1)) == 58
    assert cumulative_code((2, 2, 1, 1)) == 58
    print("FIDELITY source example: C(2,2,1,1) = 58")

    partition_codes = {}
    row_counts = []
    for total in range(1, 16):
        row = list(partitions(total))
        row_counts.append(len(row))
        assert len(row) == PARTITION_COUNTS[total - 1]
        for s in row:
            code = source_code(s)
            assert code == cumulative_code(s)
            assert decode_code(code) == s
            assert code not in partition_codes
            partition_codes[code] = s
            n = inverse_staircase(s)
            assert staircase(n) == s
            assert b(n) == code
            assert (1 << (total - 1)) <= code < (1 << total)
    assert tuple(row_counts) == PARTITION_COUNTS
    print("PARTITION_COUNTS sizes 1..15:", ",".join(map(str, row_counts)))
    print("PARTITION_CHECKS all", len(partition_codes), "partitions of sizes 1..15 passed")

    source_prefix = tuple(sorted(code for code in partition_codes if code < (1 << 8)))
    assert source_prefix == OEIS_TERMS
    print("FIDELITY complete OEIS DATA field: all", len(OEIS_TERMS),
          "terms match literal partition encoding; revision", OEIS_REVISION)
    print("OEIS_TERMS:", ",".join(map(str, source_prefix)))

    values = {}
    for n in range(1, LIMIT):
        s = staircase(n)
        assert s and all(x > 0 for x in s)
        assert all(x >= y for x, y in zip(s, s[1:]))
        assert inverse_staircase(s) == n
        assert bin(n)[2:].count("0") == n.bit_length() - n.bit_count() == s[0] - 1
        if len(s) > 1:
            remainder = n - (1 << (n.bit_length() - 1))
            assert 0 < remainder < n
            assert staircase(remainder) == s[1:]
            assert cumulative_code(s) == (1 << (s[0] - 1)) * (1 + 2 * cumulative_code(s[1:]))
        value = b(n)
        assert value == source_code(s) == cumulative_code(s)
        assert decode_code(value) == s
        assert value not in values, (n, values.get(value), value)
        values[value] = n
    print("RECURRENCE_IDENTITY, inverse staircase, code decoding, and injectivity:",
          "all 1 <= n <", LIMIT, "passed")

    expected_bounded_codes = {code for code in partition_codes if code < CODE_BOUND}
    tested_bounded_values = {value for value in values if value < CODE_BOUND}
    assert len(expected_bounded_codes) == 271
    assert tested_bounded_values == expected_bounded_codes
    assert {b(n) for n in range(1, 1 << 13) if b(n) < CODE_BOUND} == expected_bounded_codes
    print("BOUNDED_SET_EQUALITY: the tested b-values below", CODE_BOUND,
          "are exactly the", len(expected_bounded_codes), "partition codes below that bound")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
