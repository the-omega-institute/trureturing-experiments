#!/usr/bin/env python3
"""Exact finite checks for Paul D. Hanna's OEIS A382316 conjecture.

The all-index assertion is not certified by this experiment.
"""

from hashlib import sha256
from pathlib import Path


def mul(a, b):
    """Multiply integer series modulo x**len(a)."""
    assert len(a) == len(b)
    c = [0] * len(a)
    for i, ai in enumerate(a):
        if ai:
            for j in range(len(a) - i):
                c[i + j] += ai * b[j]
    return c


def subst(a, k):
    assert k >= 1
    return [a[n // k] if n % k == 0 else 0 for n in range(len(a))]


def inverse_unit(a):
    assert a[0] == 1
    b = [1] + [0] * (len(a) - 1)
    for n in range(1, len(a)):
        b[n] = -sum(a[i] * b[n - i] for i in range(1, n + 1))
    assert mul(a, b) == [1] + [0] * (len(a) - 1)
    return b


def coupled(size):
    f, g = [1] + [0] * (size - 1), [1] + [0] * (size - 1)
    for n in range(1, size):
        g[n] = int(n == 1) + sum(
            f[(n - 3 * j) // 2] * g[j]
            for j in range(n // 3 + 1) if (n - 3 * j) % 2 == 0
        )
        f[n] = g[n] + sum(f[i] * g[n - 2 * i] for i in range(1, n // 2 + 1))
    return f, g


def source_cleared(size):
    """Independent recursion from F F(x^6) = x F(x^2) F(x^6) + F(x^2)^2 F(x^3).

    This uses subtraction and never calls the positive coupled recursion.
    square[m] is the coefficient of x^m in F(x)^2.
    """
    f, square = [1] + [0] * (size - 1), [1] + [0] * (size - 1)
    for n in range(1, size):
        rhs1 = sum(f[(n - 1 - 6 * j) // 2] * f[j]
                   for j in range((n - 1) // 6 + 1)) if n % 2 else 0
        rhs2 = sum(square[(n - 3 * k) // 2] * f[k]
                   for k in range(n // 3 + 1) if (n - 3 * k) % 2 == 0)
        correction = sum(f[n - 6 * j] * f[j] for j in range(1, n // 6 + 1))
        f[n] = rhs1 + rhs2 - correction
        square[n] = sum(f[i] * f[n - i] for i in range(n + 1))
    return f


def main():
    data = Path(__file__).with_name("b382316.txt").read_bytes()
    rows = [tuple(map(int, line.split())) for line in data.decode().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    assert [n for n, _ in rows] == list(range(1, 1031))
    size = len(rows)
    f, g = coupled(size)
    assert f == [value for _, value in rows]
    assert f == source_cleared(size)
    assert all(isinstance(v, int) and v >= 1 for v in f)
    assert all(isinstance(v, int) and v >= 0 for v in g)
    assert g[0] == g[1] == 1
    assert all(f[n] >= f[n // 2] for n in range(1, size))

    f2, f3, f6, g3 = subst(f, 2), subst(f, 3), subst(f, 6), subst(g, 3)
    assert f == mul(f2, g)
    expected_g = mul(f2, g3)
    expected_g[1] += 1
    assert g == expected_g
    assert g == mul(f, inverse_unit(f2))
    assert g3 == mul(f3, inverse_unit(f6))
    rhs_unit = mul(mul(mul(f2, f2), f3), inverse_unit(f6))
    assert f == [a + b for a, b in zip([0] + f2[:-1], rhs_unit)]
    assert mul(f, f6) == [a + b for a, b in zip(
        [0] + mul(f2, f6)[:-1], mul(mul(f2, f2), f3))]

    # Literal denominator-cleared A-equation, including its initial zero coefficients.
    a = [0] + f
    a2, a3, a6 = subst(a, 2), subst(a, 3), subst(a, 6)
    assert mul(a, a6) == [u + v for u, v in zip(mul(a2, a6), mul(mul(a2, a2), a3))]

    print("Source: https://oeis.org/A382316, revision 11 (2025-05-10)")
    print("Published input: b382316.txt, sha256=" + sha256(data).hexdigest())
    print("Fidelity: all 1030 published terms a(1)..a(1030) matched exactly")
    print("First 24 terms: " + ", ".join(map(str, f[:24])))
    print("Independent scalar cleared-equation recurrence: 1030 coefficients matched")
    print("Coupled identities, quotient bridge, normalized source identity: degrees 0..1029 passed")
    print("Literal cleared A-equation: degrees 0..1030 passed")
    print("Integrality, f_n >= 1, g_n >= 0, f_n >= f_floor(n/2): degrees 0..1029 passed")
    print("Boundary: finite exact integer arithmetic; no all-index or Lean proof is claimed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
