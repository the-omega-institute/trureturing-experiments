#!/usr/bin/env python3
"""Exact finite checks for OEIS A111984 and Bala's modulo-five conjecture.

Only integer arithmetic is used. A finite check does not prove the all-index
conjecture. The two root algorithms are the cubic coefficient recurrence and
the differential identity 3*F*A' = A*F' derived from A**3 = F.
"""
from hashlib import sha256
from pathlib import Path

MAX_INDEX = 1104
DISPLAYED = (
    1, 8, -64, 960, -15360, 274432, -5223936, 103604224,
    -2115993600, 44188631040, -939004514304, 20235083907072,
    -441102457372672, 9708551114588160, -215436164616683520,
    4814308934865944576, -108243166106365722624, 2446764246113433157632,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def weighted_triangular_series(limit):
    f = [0] * (limit + 1)
    k = 0
    while k * (k + 1) // 2 <= limit:
        f[k * (k + 1) // 2] += (2 * k + 1) * 8**k
        k += 1
    return f


def cube_root(f):
    """Solve [x^n] A^3 = f[n]; maintain the square of the known prefix."""
    limit = len(f) - 1
    require(f[0] == 1, "normalization")
    a = [1] + [0] * limit
    square = [1] + [0] * limit
    for n in range(1, limit + 1):
        known_cube = square[n] + sum(a[i] * square[n - i] for i in range(1, n))
        numerator = f[n] - known_cube
        q, rem = divmod(numerator, 3)
        require(rem == 0, f"nonintegral cube-root coefficient at {n}")
        a[n] = q
        square[n] += 2 * q
        for j in range(1, min(n - 1, limit - n) + 1):
            square[n + j] += 2 * q * a[j]
        if 2 * n <= limit:
            square[2 * n] += q * q
    return a


def root_via_derivative(f):
    limit = len(f) - 1
    nonzero = [(i, x) for i, x in enumerate(f) if i and x]
    a = [1] + [0] * limit
    for n in range(1, limit + 1):
        numerator = n * f[n]
        for i, fi in nonzero:
            if i >= n:
                break
            numerator += (4 * i - 3 * n) * fi * a[n - i]
        a[n], rem = divmod(numerator, 3 * n)
        require(rem == 0, f"differential recurrence nonintegral at {n}")
    return a


def multiply(a, b, limit, modulus=None):
    c = [0] * (limit + 1)
    aa = [(i, x) for i, x in enumerate(a) if x]
    bb = [(i, x) for i, x in enumerate(b) if x]
    for i, x in aa:
        for j, y in bb:
            if i + j > limit:
                break
            c[i + j] += x * y
    if modulus is not None:
        c = [x % modulus for x in c]
    return c


def euler_product(limit):
    e = [1] + [0] * limit
    for k in range(1, limit + 1):
        for n in range(limit, k - 1, -1):
            e[n] -= e[n - k]
    return e


def main():
    limit = MAX_INDEX
    f = weighted_triangular_series(limit)
    a = cube_root(f)
    require(a[:len(DISPLAYED)] == list(DISPLAYED), "OEIS displayed terms mismatch")
    print(f"FIDELITY: matched all {len(DISPLAYED)} displayed OEIS A111984 terms (indices 0..17)")
    require(a == root_via_derivative(f), "independent root algorithms disagree")
    print(f"ROOT: independent exact recurrences agree on {len(a)} integral coefficients (0..{limit})")
    # Recompute the full defining cubic from the finished coefficient list.
    square = multiply(a, a, limit)
    require(multiply(square, a, limit) == f, "A^3 differs from F")
    print(f"DEFINING_EQUATION: A^3 = F through degree {limit}")
    e = euler_product(limit)
    ecube = multiply(multiply(e, e, limit), e, limit)
    jacobi = [0] * (limit + 1)
    k = 0
    while k * (k + 1) // 2 <= limit:
        jacobi[k * (k + 1) // 2] += (-1)**k * (2 * k + 1)
        k += 1
    require(ecube == jacobi, "truncated Jacobi identity mismatch")
    require(all((x - y) % 9 == 0 for x, y in zip(f, ecube)), "F != E^3 modulo 9")
    require(all((x - y) % 3 == 0 for x, y in zip(a, e)), "A != E modulo 3")
    print(f"INTEGRALITY_ROUTE: Jacobi, F = E^3 mod 9, A = E mod 3 through degree {limit}")
    h = [x % 5 for x in f]
    require(all(x == 0 for n, x in enumerate(h) if n % 5 not in (0, 1)), "source support")
    h2 = multiply(h, h, limit, 5)
    require(all(x == 0 for n, x in enumerate(h2) if n % 5 not in (0, 1, 2)), "squared support")
    b = [x % 5 for x in a]
    bx5 = [b[n // 5] if n % 5 == 0 else 0 for n in range(limit + 1)]
    require(multiply(b, bx5, limit, 5) == h2, "B*B(x^5) != H^2")
    print(f"SUPPORT: H in classes 0,1; H^2 in 0,1,2; B*B(x^5) = H^2 through degree {limit}")
    for residue in (4, 3):
        indices = list(range(residue, limit + 1, 5))
        require(all(a[n] % 5 == 0 for n in indices), f"failed residue {residue}")
        print(f"CONGRUENCE: 5 divides a(5*n+{residue}) for n=0..{(indices[-1]-residue)//5} ({len(indices)} cases)")
    require(a[0] % 5 != 0 and a[1] % 5 != 0 and a[2] % 5 != 0, "negative controls")
    print("CONTROLS: each proposed extra zero class 0,1,2 fails already at its first index")
    digest = sha256(("\n".join(map(str, a)) + "\n").encode()).hexdigest()
    print(f"COEFFICIENT_SHA256 (one decimal integer per LF-terminated line): {digest}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
