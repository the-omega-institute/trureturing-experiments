#!/usr/bin/env python3
"""Exact checks for OEIS A340991's shifted-array recursion."""

from __future__ import annotations


def poly_add(a: list[int], b: list[int]) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i] += v
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_mul(a: list[int], b: list[int]) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def poly_scale(a: list[int], c: int) -> list[int]:
    return [c * x for x in a]


def poly_shift(a: list[int]) -> list[int]:
    return [0] + a


def primes(count: int) -> list[int]:
    out: list[int] = []
    n = 2
    while len(out) < count:
        d = 2
        prime = True
        while d * d <= n:
            if n % d == 0:
                prime = False
                break
            d += 1
        if prime:
            out.append(n)
        n += 1
    return out


def row_from_coefficients(a: list[int], n: int) -> list[int]:
    """T(n,r) from the finite Cauchy product, for 0 <= r <= n."""
    if n == 0:
        return [1]
    row = [0] * (n + 1)
    for r in range(1, n + 1):
        # DP for compositions of n into r positive parts.
        dp = {(0, 0): 1}
        for _ in range(r):
            nxt: dict[tuple[int, int], int] = {}
            for (parts, total), value in dp.items():
                for part in range(1, n - total + 1):
                    key = (parts + 1, total + part)
                    nxt[key] = nxt.get(key, 0) + value * a[part]
            dp = nxt
        row[r] = dp.get((r, n), 0)
    return row


def kurkov_rows(a: list[int], max_n: int, max_k: int) -> dict[tuple[int, int], list[int]]:
    """Compute R(n,k) as exact integer polynomials using the stated recursion."""
    rows: dict[tuple[int, int], list[int]] = {}
    for k in range(1, max_k + max_n + 2):
        rows[(1, k)] = [a[k]]
    for n in range(2, max_n + 1):
        for k in range(1, max_k + max_n - n + 2):
            first = rows[(n - 1, k + 1)]
            second = poly_shift(poly_scale(rows[(n - 1, 1)], a[k]))
            rows[(n, k)] = poly_add(first, second)
    return rows


def check_sequence(a: list[int], max_n: int = 10) -> None:
    rows = kurkov_rows(a, max_n, max_n)
    for n in range(1, max_n + 1):
        coeff_row = row_from_coefficients(a, n)
        assert poly_shift(rows[(n, 1)]) == coeff_row, (n, rows[(n, 1)], coeff_row)
        for k in range(1, min(5, max_n) + 1):
            # Independent composition expansion for Q(n,k), including the
            # shifted final factor a_{c_l+k-1}.
            q = [0]
            for length in range(1, n + 1):
                q = poly_add(q, composition_weight(a, n, length, k))
            assert rows[(n, k)] == q, (n, k, rows[(n, k)], q)


def composition_weight(a: list[int], total: int, length: int, k: int) -> list[int]:
    """Q(total,k): sum x^(length-1) product before shifted last factor."""
    total_weight = 0

    def visit(remaining: int, used: int, product: int) -> None:
        nonlocal total_weight
        if used == length - 1:
            final_part = remaining
            if final_part >= 1:
                total_weight += product * a[final_part + k - 1]
            return
        for part in range(1, remaining - (length - used - 1) + 1):
            visit(remaining - part, used + 1, product * a[part])

    visit(total, 0, 1)
    return [0] * (length - 1) + [total_weight]


def main() -> None:
    ps = [0] + primes(40)
    check_sequence(ps, 10)
    # A signed and repeated sequence checks the claimed arbitrary-input
    # extension without relying on primality or positivity.
    arbitrary = [0] + [((7 * i + 3) % 11) - 5 for i in range(1, 45)]
    check_sequence(arbitrary, 9)
    row4 = row_from_coefficients(ps, 4)
    assert row4 == [0, 7, 29, 36, 16], row4
    print("OEIS A340991 row 4: [0, 7, 29, 36, 16]")
    print("Checked universal recursion for primes and a signed arbitrary sequence through n=10.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
