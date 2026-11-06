#!/usr/bin/env python3
"""Exact checks for the OEIS A359506 sharp 5/3 endpoint conjecture."""


def endpoint_bruteforce(n: int, limit: int = 10000) -> int:
    """Return the least endpoint by the OEIS incremental XOR-reachability algorithm."""
    if n == 0:
        return 0
    reachable = {n}
    for m in range(n + 1, limit + 1):
        old = tuple(reachable)
        if m in reachable:
            return m
        reachable.update(x ^ m for x in old)
    raise AssertionError(f"no endpoint found below {limit} for n={n}")


def lowbit(x: int) -> int:
    return x & -x


def exact_formula(n: int) -> int:
    if n == 0:
        return 0
    q = lowbit(n + 1)
    u = (n + 1) // q
    if u == 1:
        return 3 * q // 2
    return q * (u + 2)


def main() -> None:
    # Fidelity: the first examples in OEIS A359506 (including the zero case).
    expected_prefix = [0, 3, 5, 6, 7, 10, 9, 12, 11, 14, 13]
    observed_prefix = [endpoint_bruteforce(n) for n in range(len(expected_prefix))]
    assert observed_prefix == expected_prefix, (observed_prefix, expected_prefix)
    print("fidelity prefix A359506[0..10] =", observed_prefix)

    # Independent exact minimum checks on every n in a nontrivial small range.
    for n in range(0, 513):
        brute = endpoint_bruteforce(n)
        formula = exact_formula(n)
        assert brute == formula, (n, brute, formula)
    print("exact minima checked for n = 0..512")

    # The claimed bound, using integer arithmetic only.
    for n in range(0, 513):
        e = exact_formula(n)
        assert 3 * e <= 5 * (n + 1), (n, e)
    print("3*E(n) <= 5*(n+1) checked for n = 0..512")

    # Equality family n = 3*2^s - 1, E(n) = 5*2^s.
    for s in range(0, 13):
        n = 3 * (1 << s) - 1
        e = exact_formula(n)
        assert e == 5 * (1 << s), (s, n, e)
        assert 3 * e == 5 * (n + 1), (s, n, e)
    print("sharp equality family checked for s = 0..12")

    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
