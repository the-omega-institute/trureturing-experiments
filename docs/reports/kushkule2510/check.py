"""Exact checks for arXiv:2510.04458v1."""

from collections import Counter
from fractions import Fraction
from functools import cache
from itertools import combinations, product


Element = tuple[int, int]
Trace = tuple[int, ...]


def multiply(m: int, x: Element, y: Element) -> Element:
    """Multiply r^j s^e in coordinates (j, e)."""
    j, e = x
    k, f = y
    return ((j + (-1)**e * k) % m, (e + f) % 2)


def inverse(m: int, x: Element) -> Element:
    j, e = x
    return ((-(-1)**e * j) % m, e)


def commutator(m: int, x: Element, y: Element) -> Element:
    return multiply(m, multiply(m, multiply(m, x, y), inverse(m, x)), inverse(m, y))


def trace(k: int, m: int, g: Element) -> Trace:
    """Represent the real trace by the exponents of its root-of-unity terms."""
    j, e = g
    if k == 0:
        return (0,)  # Trivial character.
    if k == -1:
        return (0,)  # The sign character is handled by its coefficient.
    if e:
        return ()
    return ((k * j) % m, (-k * j) % m)


def sign(k: int, g: Element) -> int:
    return -1 if k == -1 and g[1] else 1


def divide_monic(
    coefficients: list[Fraction], divisor: tuple[Fraction, ...]
) -> tuple[list[Fraction], list[Fraction]]:
    """Divide polynomials with coefficients in ascending degree order."""
    remainder = coefficients.copy()
    degree = len(divisor) - 1
    quotient = [Fraction(0)] * max(0, len(remainder) - degree)
    for shift in reversed(range(len(quotient))):
        leading = remainder[degree + shift]
        quotient[shift] = leading
        for i, coefficient in enumerate(divisor):
            remainder[i + shift] -= leading * coefficient
    return quotient, remainder[:degree]


@cache
def cyclotomic(m: int) -> tuple[Fraction, ...]:
    """Use X^m - 1 = ∏_{d | m} Φ_d(X)."""
    coefficients = [Fraction(-1)] + [Fraction(0)] * (m - 1) + [Fraction(1)]
    for divisor in range(1, m):
        if m % divisor == 0:
            coefficients, remainder = divide_monic(coefficients, cyclotomic(divisor))
            assert not any(remainder)
    return tuple(coefficients)


def evaluate_equal(m: int, terms: list[tuple[Fraction, Trace]]) -> bool:
    """Certify a polynomial vanishes at a primitive m-th root."""
    coefficients = [Fraction(0)] * m
    for coefficient, exponents in terms:
        for exponent in exponents:
            coefficients[exponent % m] += coefficient
    _, remainder = divide_monic(coefficients, cyclotomic(m))
    return not any(remainder)


def check_quaternions() -> None:
    # Basis 1, i, j, k. Signs are folded into each table entry.
    table = ((1, 2, 3, 4), (2, -1, 4, -3),
             (3, -4, -1, 2), (4, 3, -2, -1))
    elements = [unit * sign_value for unit in range(1, 5) for sign_value in (1, -1)]

    def quaternion_mul(x: int, y: int) -> int:
        signed_entry = table[abs(x) - 1][abs(y) - 1]
        return signed_entry if (x > 0) == (y > 0) else -signed_entry

    count = sum(quaternion_mul(x, y) == quaternion_mul(y, x)
                for x, y in product(elements, repeat=2))
    probability = Fraction(count, len(elements) ** 2)
    assert probability == Fraction(5, 8)
    print(f"Source Example 1, Q8: {count}/64 commuting pairs = {probability}")


def check_dihedral(m: int) -> None:
    elements = list(product(range(m), (0, 1)))
    identity = (0, 0)
    order = len(elements)
    counts = Counter(commutator(m, x, y) for x, y in product(elements, repeat=2))
    assert sum(counts.values()) == order**2

    irreducibles = {-1: 1, 0: 1, **{k: 2 for k in range(1, (m + 1) // 2)}}
    assert sum(d**2 for d in irreducibles.values()) == order
    for i, j in combinations(irreducibles, 2):
        inner_product = [
            (Fraction(sign(i, g) * sign(j, g)),
             tuple(x - y for x, y in product(trace(i, m, g), trace(j, m, g))))
            for g in elements
        ]
        assert evaluate_equal(m, inner_product)
    for i in irreducibles:
        squared_norm = [
            (Fraction(1), tuple(x - y for x, y in product(trace(i, m, g), repeat=2)))
            for g in elements
        ]
        assert evaluate_equal(m, squared_norm + [(Fraction(-order), (0,))])

    probability = {g: Fraction(counts[g], order**2) for g in elements}
    assert probability[identity] == Fraction(m + 3, 4 * m)
    assert probability[identity] == Fraction(len(irreducibles), order)  # (1.2)
    lower_gaps = []
    for g in elements:
        exponent, reflection = g
        if reflection:
            assert probability[g] == 0
            lower_rhs = Fraction(1, 3)
        elif exponent:
            assert probability[g] == Fraction(3, 4 * m)
            lower_rhs = Fraction(1, 6)  # 1 + (ζ^j + ζ^-j)/4 >= 1/2.
        else:
            lower_rhs = Fraction(1, 2)
        lower_gaps.append(lower_rhs - probability[g])

        frobenius = [
            (Fraction(sign(k, g), order * degree), trace(k, m, g))
            for k, degree in irreducibles.items()
        ]
        assert evaluate_equal(m, frobenius + [(-probability[g], (0,))])  # (1.1)
        regular = [
            (Fraction(degree * sign(k, g)), trace(k, m, g))
            for k, degree in irreducibles.items()
        ]
        assert evaluate_equal(m, regular + [(Fraction(-order * (g == identity)), (0,))])

        for b in (0, 1, 3):
            coefficients = {k: (k in (0, 1)) + b * d for k, d in irreducibles.items()}
            if b:
                assert all(a > 0 for a in coefficients.values())
            weighted_trace = [
                (Fraction(coefficients[k] * sign(k, g), d**2), trace(k, m, g))
                for k, d in irreducibles.items()
            ]
            averaged_trace = []
            for h, count in counts.items():
                t = multiply(m, inverse(m, h), g)
                averaged_trace.extend(
                    (Fraction(count * a * sign(k, t), order**2), trace(k, m, t))
                    for k, a in coefficients.items()
                )
            assert evaluate_equal(m, averaged_trace + [(-c, e) for c, e in weighted_trace])

            virtual_degree = sum(coefficients[k] * d for k, d in irreducibles.items())
            assert virtual_degree == 3 + order * b
            weighted_gap = weighted_trace + [(-virtual_degree * probability[g], (0,))]
            base_gap = [(Fraction(1) - 3 * probability[g], (0,)),
                        (Fraction(1, 4), trace(1, m, g))]
            assert evaluate_equal(m, weighted_gap + [(-c, e) for c, e in base_gap])

    witness = ((m - 1) // 2, 0)
    assert witness != identity
    assert trace(1, m, witness) == ((m - 1) // 2, (m + 1) // 2)
    # φ(witness) = 1 - 2 cos(π/m). Use cos x >= 1 - x²/2 and π < 22/7.
    upper_witness = -1 + Fraction(22, 7 * m) ** 2
    assert upper_witness < 0
    assert min(lower_gaps) > 0
    print(f"m={m:2}: c(1)={probability[identity]}, "
          f"minimum gap >= {min(lower_gaps)}, witness <= {upper_witness}")


def main() -> None:
    check_quaternions()
    for m in range(5, 22, 2):
        check_dihedral(m)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
