"""Exact checks for Torun, arXiv:2611.16719v1, Proposition 2."""
from fractions import Fraction
from itertools import combinations
from random import Random

from sympy import Poly, Rational, cancel, expand, factor, symbols


N = 4
SIZE = 1 << N
PARITY = SIZE - 1
POINTS = range(SIZE)
WEIGHTS = tuple(k.bit_count() for k in POINTS)


def character(k, x):
    return (-1) ** ((k & x).bit_count())


WALSH = tuple(tuple(character(k, x) for x in POINTS) for k in POINTS)


def correlations(signs):
    return tuple(sum(a * b for a, b in zip(row, signs, strict=True)) for row in WALSH)


def check_source():
    """Verify the exact branch and its matching Proposition 2 upper bound."""
    r = symbols("r", positive=True)
    s = (1 - r) / (1 + r)
    for sigma in (1, -1):
        diag = 1 / (2 * (1 + sigma * s))
        off_diag = sigma * s * diag
        assert cancel(diag + off_diag) == Rational(1, 2)
        assert cancel(diag - off_diag) == r**sigma / 2

    # Appendix C matrices, encoded by their columns (low bit is row one).
    for n, columns in ((2, (0b01, 0b11)), (3, (0b011, 0b101, 0b111))):
        images = []
        for k in range(1 << n):
            image = 0
            for i, col in enumerate(columns):
                if k & (1 << i):
                    image ^= col
            images.append(image)
        assert len(set(images)) == 1 << n
        exponents = [k.bit_count() + img.bit_count() for k, img in enumerate(images)]
        assert max(exponents) == n + 1
        # Top source character: any permuted image is balanced, so the target
        # quadratic form is bounded below by the smallest nonconstant eigenvalue.
        assert min(k.bit_count() for k in range(1, 1 << n)) == 1
        assert cancel((r**n / 2**n) / (r**-1 / 2**n)) == r ** (n + 1)
        print(f"Proposition 2: N={n}, exact optimum r^{n + 1}")


def check_balanced():
    for k in POINTS:
        assert correlations(WALSH[k]) == tuple(SIZE if j == k else 0 for j in POINTS)
    counted = dictators = 0
    minimum = SIZE**2
    for subset in combinations(POINTS, SIZE // 2):
        positive = frozenset(subset)
        signs = tuple(1 if x in positive else -1 for x in POINTS)
        values = correlations(signs)
        assert values[0] == 0
        assert all(value % 4 == 0 for value in values)
        assert sum(value**2 for value in values) == SIZE**2
        high_degree = sum(value**2 for k, value in enumerate(values) if WEIGHTS[k] >= 2)
        if high_degree == 0:
            dictators += 1
            assert sum(abs(values[k]) == SIZE for k in POINTS if WEIGHTS[k] == 1) == 1
        else:
            minimum = min(minimum, high_degree)
            assert high_degree >= 16
        counted += 1
    assert (counted, dictators, minimum) == (12870, 8, 64)
    print(f"Balanced signs: {counted}, signed coordinates: {dictators}, "
          f"minimum positive high-degree mass: {Fraction(minimum, SIZE**2)}")


def certificate_coefficients(r):
    a = 256 * r**2 / (1 + 15 * r)
    b = 16 * r**6 * (1 - r) / (1 + 15 * r)
    return a, b


def check_family():
    r, u, v = symbols("r u v", positive=True)
    a, b = certificate_coefficients(r)
    min_mass = Rational(1, 16)
    first_case = a * (r**-1 + min_mass * (r**-2 - r**-1)) / SIZE
    assert cancel(first_case) == 1
    for j in (1, 2, 3):
        # u = 1/r = 1 + v. A polynomial with nonnegative coefficients
        # at v >= 0 verifies each retained spectral weight.
        direct = u ** (4 + j) if 2 <= j <= 3 else 0
        paired = u ** (3 + (4 - j)) if 2 <= 4 - j <= 3 else 0
        assert all(
            coeff >= 0
            for coeff in Poly(expand((direct + paired - u**6).subs(u, 1 + v)), v).all_coeffs()
        )
    assert cancel((a / r + b / r**6) / SIZE) == 1
    objective = (a * r**4 + sum(b for j in WEIGHTS if 2 <= j <= 3)) / SIZE
    assert cancel(objective - r**6 * (26 - 10 * r) / (1 + 15 * r)) == 0
    gap = factor(r**5 - objective)
    assert cancel(gap - r**5 * (1 - r) * (1 - 10 * r) / (1 + 15 * r)) == 0
    # After r = 1/(10+v), the gap numerator has positive coefficients.
    numerator, denominator = cancel(gap.subs(r, 1 / (10 + v))).as_numer_denom()
    for expr in (numerator, denominator):
        coefficients = Poly(expr, v).all_coeffs()
        assert any(c > 0 for c in coefficients) and all(c >= 0 for c in coefficients)
    print("Both dual-case lower bounds equal 1. The gap is positive for 0 < r < 1/10.")


def dual_trace(permutation, dual_eigenvalues, target_eigenvalues):
    value = Fraction(0)
    for k, eigenvalue in dual_eigenvalues.items():
        transformed = tuple(WALSH[k][permutation[x]] for x in POINTS)
        squares = tuple(Fraction(c, SIZE)**2 for c in correlations(transformed))
        value += eigenvalue * sum(y * d for y, d in zip(squares, target_eigenvalues, strict=True))
    return value


def check_rational_instance():
    r = Fraction(1, 16)
    s = Fraction(15, 17)
    assert Fraction(8, 17)**2 + s**2 == 1 and r == (1 - s) / (1 + s)
    a, b = certificate_coefficients(r)
    dual = {PARITY: a} | {k: b / r**w for k, w in enumerate(WEIGHTS) if 2 <= w <= 3}
    assert a == Fraction(16, 31)
    assert b / r**2 == Fraction(15, 126976)
    assert b / r**3 == Fraction(15, 7936)
    source = tuple(r**w / SIZE for w in WEIGHTS)
    target = tuple(r**-w / SIZE for w in WEIGHTS)
    for x in POINTS:
        for y in POINTS:
            distance = (x ^ y).bit_count()
            assert sum(source[k] * WALSH[k][x] * WALSH[k][y] for k in POINTS) / SIZE == (
                Fraction(17, 64) ** (4 - distance) * Fraction(15, 64) ** distance
            )
            assert sum(target[k] * WALSH[k][x] * WALSH[k][y] for k in POINTS) / SIZE == (
                Fraction(17, 4) ** (4 - distance) * Fraction(-15, 4) ** distance
            )

    random = Random(261116719)
    for _ in range(64):
        permutation = random.sample(POINTS, SIZE)
        assert dual_trace(permutation, dual, target) >= 1

    # Use arbitrary bijections within the two sign fibers to force Case 2.
    for i in (1, 2, 4, 8):
        for sign in (1, -1):
            for _ in range(8):
                permutation = [None] * SIZE
                for value in (1, -1):
                    domain = [x for x in POINTS if WALSH[i][x] == sign * value]
                    image = [x for x in POINTS if WALSH[PARITY][x] == value]
                    random.shuffle(image)
                    for x, y in zip(domain, image, strict=True):
                        permutation[x] = y
                assert len(set(permutation)) == SIZE
                assert all(WALSH[PARITY][permutation[x]] == sign * WALSH[i][x] for x in POINTS)
                assert dual_trace(permutation, dual, target) >= 1

    obj = sum(dual[k] * source[k] for k in dual)
    assert obj == Fraction(203, 260046848) == Fraction(203, 248) * r**5 < r**5
    print("Four-copy dual trace: 64 sampled bijections and 64 signed-coordinate bijections.")
    print(f"r={r}, s={s}: dual objective {obj} < {r**5}")


if __name__ == "__main__":
    check_source()
    check_balanced()
    check_family()
    check_rational_instance()
    print("ALL_CHECKS_OK")
