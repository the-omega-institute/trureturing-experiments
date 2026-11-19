"""Exact coefficient checks for the conjecture after Chen–Mortici (50)."""

from fractions import Fraction
from math import comb, factorial


def bernoulli_numbers(last):
    """Use (11) after multiplying both sides by exp(z) - 1."""
    numbers = [Fraction(1)]
    for n in range(1, last + 1):
        previous = sum(comb(n + 1, k) * numbers[k] for k in range(n))
        numbers.append(-previous / (n + 1))
    return numbers


def coefficient(j, bernoulli):
    """Compute the coefficient in (45) for j >= 1."""
    elementary = Fraction((-1) ** j, j + 1)
    if j == 1:
        return elementary
    return elementary + bernoulli[j] / (j * (j - 1))


def exponential_coefficients(exponent):
    """Use (exp P)' = P' exp P for the formal series with P(0) = 0."""
    result = [Fraction(1)]
    for n in range(1, len(exponent)):
        result.append(sum(j * exponent[j] * result[n - j]
                          for j in range(1, n + 1)) / n)
    return result


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    bernoulli = bernoulli_numbers(40)
    omega = [Fraction(0)] + [coefficient(j, bernoulli) for j in range(1, 41)]

    expected_log = [Fraction(-1, 2), Fraction(5, 12), Fraction(-1, 4)]
    expected_exp = [Fraction(1), Fraction(-1, 2), Fraction(13, 24), Fraction(-23, 48)]
    require(omega[1:4] == expected_log, "Equation (49) differs")
    require(exponential_coefficients(omega[:4]) == expected_exp, "Equation (50) differs")
    print("Equation (49):", ", ".join(map(str, expected_log)))
    print("Equation (50):", ", ".join(map(str, expected_exp)))

    # pi = 4 arctan(1). An even-length alternating sum is a strict
    # lower bound. Its one-term truncation gives the strict upper bound 4.
    pi_lower = 4 * sum(Fraction((-1) ** k, 2 * k + 1) for k in range(8))
    require(pi_lower > 3, "Failed to certify pi > 3")

    # Gamma(2) = 1!, so at x = 1 the squared log argument is 2/pi.
    # log(u) < u - 1 for 0 < u < 1 gives the strict bound used here.
    squared_numerator = (2 * factorial(1)) ** 2
    log_argument_square_upper = Fraction(squared_numerator, 2 * 3)
    h_upper = (log_argument_square_upper - 1) / 2
    require(h_upper == Fraction(-1, 6), "Wrong bound on H(1)")
    partial_sum = sum(omega[1:20])
    require(partial_sum == Fraction(-663783049, 6983776800), "Wrong partial sum")
    require(partial_sum > Fraction(-1, 10), "Partial sum insufficient")
    require((-1) ** (19 - 1) == 1, "Wrong sign of G_19")
    require(h_upper - partial_sum < Fraction(-1, 15), "Counterexample bound failed")
    print("sum(omega_1, ..., omega_19) =", partial_sum)
    print("Certified: G_19(1) < -1/15")

    # Euler's even-zeta formula gives the ratio A_(n+4)/A_n.
    # For n >= 20 each of the four factorial-ratio factors is increasing.
    # Use pi < 4 and 1 < zeta(n+4), zeta(n) < 2 in the ratio bound.
    a20 = -bernoulli[20] / (20 * 19)
    require(a20 == Fraction(174611, 125400), "Wrong A_20")
    require(a20 > Fraction(1, 21), "Failed base bound")
    ratio_lower = Fraction(19 * 20 * 21 * 22, 2 * 8 ** 4)
    require(ratio_lower > 1, "Failed uniform induction factor")
    for n in range(20, 41, 4):
        require(omega[n] < 0, f"Wrong sign for omega_{n}")
        require((-1) ** (n - 2) == 1 and (-1) ** (n - 1) == -1,
                f"Wrong limit sign for the orders {n - 1}, {n}")
    print("A_20 =", a20)
    print("Uniform lower bound on the A_(n+4)/A_n ratio for n >= 20:", ratio_lower)
    print("Exact coefficient signs for n = 20, 24, 28, 32, 36, 40: strictly negative")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
