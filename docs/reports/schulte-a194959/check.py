#!/usr/bin/env python3
"""Exact checks for Schulte's last FORMULA conjecture in OEIS A194959.

Source: https://oeis.org/A194959, revision 53; weight: A209229, revision 116.
Published data/formulas are credited to OEIS and their named contributors.
This finite check is evidence of fidelity, not a proof for every positive n.
"""

SOURCE_PREFIX = [1, 1, 2, 1, 3, 2, 1, 3, 4, 2, 1, 3, 5, 4, 2, 1, 3, 5, 6, 4, 2, 1, 3, 5, 7, 6, 4, 2, 1, 3, 5, 7, 8, 6, 4, 2, 1, 3, 5, 7, 9, 8, 6, 4, 2, 1, 3, 5, 7, 9, 10, 8, 6, 4, 2, 1, 3, 5, 7, 9, 11, 10, 8, 6, 4, 2, 1, 3, 5, 7, 9, 11, 12, 10, 8, 6, 4, 2, 1, 3, 5, 7, 9, 11, 13, 12, 10, 8, 6, 4, 2, 1, 3, 5]
LIMIT = 4096


def require(ok, label):
    if not ok:
        raise AssertionError(label)


def array_entry(r, c):
    require(r >= 1 and c >= 1, "positive array indices")
    return 2 * c - 1 if c <= r else 2 * r


def weight(k):
    return int(k > 0 and (k & (k - 1)) == 0)


def main():
    powers = set()
    power = 1
    while power <= LIMIT:
        powers.add(power)
        power *= 2
    for k in range(LIMIT + 1):
        require(weight(k) == int(k in powers), ("weight", k))
    require(weight(0) == 0 and weight(1) == 1, "A209229 offset and 2^0")

    row = []
    prefix = []
    entry_comparisons = 0
    reversed_values = []
    orientation_control = None
    for n in range(1, LIMIT + 1):
        # One-based insertion position p(n)=1+floor(n/2); Python is zero-based.
        row.insert(n // 2, n)
        if len(prefix) < len(SOURCE_PREFIX):
            prefix.extend(row[:len(SOURCE_PREFIX) - len(prefix)])
        for c, entry in enumerate(row, 1):
            r = n + 1 - c
            require(entry == array_entry(r, c), ("source array formula", n, c))
            entry_comparisons += 1
        require(sum(row) == n * (n + 1) // 2, ("published row-sum fidelity", n))

        # A(k,n+1-k) is position n+1-k of row g(n), hence row[n-k].
        actual = sum(row[n - k] * weight(k) for k in range(1, n + 1))
        from_formula = sum(array_entry(k, n + 1 - k) * weight(k)
                           for k in range(1, n + 1))
        require(actual == from_formula == 2 * n - 1, ("target", n, actual, from_formula))
        if n <= 12:
            reversed_values.append(actual)
        if n == 3:
            orientation_control = sum(row[k - 1] * weight(k) for k in range(1, n + 1))
            require(orientation_control == 4 and actual == 5, "reversal sensitivity")

    require(prefix == SOURCE_PREFIX, "published A194959 prefix fidelity")
    print("Arithmetic: Python integers only")
    print(f"Published A194959 prefix: {len(prefix)} terms matched")
    print(f"A209229 indicator: k=0..{LIMIT} matched explicit powers; weight(0)=0, weight(1)=1")
    print(f"Insertion rows vs published square-array formula: {entry_comparisons} entries matched")
    print(f"Published unweighted row-sum identity: n=1..{LIMIT}, zero failures")
    print(f"Reversed-antidiagonal conjecture: n=1..{LIMIT}, zero failures")
    print("Target values n=1..12: " + ",".join(map(str, reversed_values)))
    print(f"Orientation control n=3: forward={orientation_control}, required reversed=5")
    print("Finite verification only; no universal proof or Lean-kernel settlement is claimed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
