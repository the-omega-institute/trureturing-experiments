"""Exact endpoint signs for a four-layer ambiguity in the scalar portal response."""
from fractions import Fraction
from math import lcm, prod


def logarithmic_sign(masses, differences):
    """Sign of sum k*m²*log(m), assuming the weighted sum of m² is zero."""
    common = lcm(*(mass.denominator for mass in masses))
    integers = [int(mass * common) for mass in masses]
    powers = [digit * n * n for digit, n in zip(differences, integers, strict=True)]
    assert sum(powers) == 0
    numerator = prod(n ** power for n, power in zip(integers, powers, strict=True) if power > 0)
    denominator = prod(n ** (-power) for n, power in zip(integers, powers, strict=True) if power < 0)
    return (numerator > denominator) - (numerator < denominator)


def main():
    # At normalized slope 1, digits 1 and 3 select the uncolored and colored
    # species respectively. All remaining species take normalized slope 2.
    first = (3, 1, 3, 0)
    second = (1, 3, 1, 1)
    differences = tuple(a - b for a, b in zip(first, second, strict=True))
    for x, expected in ((Fraction(1, 10), -1), (Fraction(1, 8), 1)):
        half = Fraction(1, 2)
        gap = half - x
        pair_sum = (half - x*x) / gap
        masses = (x, (pair_sum - gap)/2, (pair_sum + gap)/2, Fraction(1))
        assert 0 < masses[0] < masses[1] < masses[2] < masses[3]
        for power in (1, 2):
            assert sum(k*m**power for k, m in zip(differences, masses, strict=True)) == 0
        sign = logarithmic_sign(masses, differences)
        assert sign == expected
        print("squared masses:", *(str(m) for m in masses), "logarithmic sign:", sign)


if __name__ == "__main__":
    main()
