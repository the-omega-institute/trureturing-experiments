#!/usr/bin/env python3
"""Exact checks for the A103528 floor-sum difference conjecture."""


def a(n: int) -> int:
    return sum(2 ** (k - 1) for k in range(1, n) if n % (2**k) == k)


def floor_log2(n: int) -> int:
    if n < 1:
        raise ValueError("floor_log2 is used only for positive n")
    return n.bit_length() - 1


def c(n: int, m: int) -> int:
    return n - n % (2**m)


def b(n: int) -> int:
    return sum(c(n - k, k) for k in range(floor_log2(n) + 1))


def main() -> None:
    expected_prefix = [0, 0, 1, 0, 1, 2, 1, 0, 1, 2, 5]
    observed_prefix = [a(n) for n in range(1, 12)]
    assert observed_prefix == expected_prefix
    assert a(1) == 0

    checked = 0
    for n in range(2, 4097):
        lhs = 2 * a(n)
        rhs = b(n) - b(n - 1) - 1
        assert lhs == rhs, (n, lhs, rhs, b(n), b(n - 1))
        checked += 1

    print("fidelity A103528 n=1..11:", observed_prefix)
    print("identity cases checked:", checked)
    for n in (2, 11, 20, 4096):
        print(f"n={n}: A={a(n)}, B={b(n)}, B_previous={b(n - 1)}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
