#!/usr/bin/env python3
"""Exact integer check for OEIS A347205's finite row conjecture."""


def v2(n: int) -> int:
    if n <= 0:
        return 0
    r = 0
    while n % 2 == 0:
        n //= 2
        r += 1
    return r


def scalar_values(limit: int) -> list[int]:
    a = [0] * (limit + 1)
    a[0] = 1
    for n in range(1, limit + 1):
        if n & 1:
            a[n] = a[n // 2]
        else:
            q = n // 2
            a[n] = a[q] + a[q - (1 << v2(q))]
    return a


def rows(limit: int) -> list[list[int]]:
    row = [[1]]
    for n in range(1, limit + 1):
        prev = row[n // 2]
        if n & 1:
            row.append(prev + [prev[-1]])
        else:
            acc = 0
            out = []
            for value in prev:
                acc += value
                out.append(acc)
            row.append(out)
    return row


def popcount(n: int) -> int:
    return n.bit_count()


def main() -> None:
    # The first terms reproduce the OEIS data field exactly through n = 15.
    expected_prefix = [1, 1, 2, 1, 3, 2, 3, 1, 4, 3, 5, 2, 6, 3, 4, 1]
    limit = 16_384
    a = scalar_values(limit * 2)
    r = rows(limit)
    assert a[:16] == expected_prefix, (a[:16], expected_prefix)

    for n in range(limit + 1):
        row = r[n]
        assert len(row) == popcount(n) + 1, (n, len(row), popcount(n))
        assert a[n] == row[-1], (n, a[n], row[-1])
        assert a[2 * n] == sum(row), (n, a[2 * n], sum(row))

    # Check the scalar recurrence independently on every positive even index.
    for n in range(1, limit + 1, 2):
        m = 2 * n
        assert a[m] == a[n] + a[n - (1 << v2(n))]

    print(f"checked scalar and row identities for n=0..{limit}")
    print(f"reproduced OEIS prefix: {expected_prefix}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
