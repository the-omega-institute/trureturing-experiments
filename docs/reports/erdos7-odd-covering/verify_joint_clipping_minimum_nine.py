#!/usr/bin/env python3
"""Replay the joint clipping bound with moduli 15 and 21 absent."""

from fractions import Fraction
from verify_minimum_modulus_nine import BANDS as CLIPPING, T as END
from verify_rh_prime_continuation import SCALE, ceiling as ceil, prime_sieve, verify_quadratic

WIDTH = 1850


def upper_weights(prime):
    tail = [0, ceil(Fraction(SCALE, prime))]
    point = [0, 0]
    power = prime
    for j in range(2, WIDTH + 1):
        power *= prime
        tail.append(ceil(Fraction(SCALE, power)))
        point.append(ceil(Fraction(SCALE * (prime - 1), power)))
        if point[-1] == 1:
            missing = WIDTH - j
            tail.extend([1] * missing)
            point.extend([1] * missing)
            break
    return tail, point


def positive_multipliers(value, prime, upper, polynomial, weights):
    tail, point = weights
    m = WIDTH // value
    a, b, c = polynomial
    geometric_mean = Fraction(prime, prime - 1)
    geometric_second = Fraction(prime * (prime + 1), (prime - 1)**2)
    high = (a * value**2 * (m*m + 2*m*geometric_mean + geometric_second)
            + b * value * (m + geometric_mean) + c)
    if high < 0:
        raise ValueError("negative upper bound for the high-state contribution")
    low = sum(point[i] * upper[i*value] for i in range(2, m + 1))
    return ceil((low + tail[m] * high) / SCALE)


def expectation(at_one, positive, cap, prime):
    return (1 - cap / prime) * at_one + cap * positive


def extend(prime, clipping, exclusion, upper, polynomial):
    pure_budget = Fraction(1, prime) if prime in (3, 5, 7) else Fraction(1)
    d = prime - 1 - pure_budget
    rho = (prime - 1) / d
    g = 1 / (1 - clipping)
    cap = rho * g
    threshold = 1 + exclusion + clipping * d
    if not (threshold <= WIDTH and cap <= prime):
        raise ValueError(f"invalid width or cap at {prime}")
    weights = upper_weights(prime)
    output = [0]
    for value in range(1, WIDTH + 1):
        positive = positive_multipliers(value, prime, upper, polynomial, weights)
        ordinary = expectation(upper[value], positive, cap, prime)
        if value < threshold:
            actual_fraction_bound = (value - 1 - exclusion) / d
            conditional_cap = rho / (1 - max(0, actual_fraction_bound))
            ordinary = max(
                expectation(upper[value], positive, conditional_cap, prime),
                ordinary - SCALE * (cap - conditional_cap) / cap,
            )
        loss = max(0, value - threshold) * g / d
        output.append(ceil(ordinary + SCALE * loss))
    a, b, c = polynomial
    high_polynomial = (
        ceil(a * (1 + cap*Fraction(3*prime - 1, (prime - 1)**2))),
        ceil(b * (1 + cap/(prime - 1)) + SCALE*g/d),
        ceil(c - SCALE*g*threshold/d),
    )
    return output, high_polynomial


def run():
    sieve = prime_sieve(5_000_000)
    tail_bound = verify_quadratic(sieve)
    upper = [tail_bound * value**2 for value in range(WIDTH + 1)]
    polynomial = tail_bound, 0, 0
    count = 0
    for prime in range(END, 7, -1):
        if not sieve[prime]:
            continue
        clipping = Fraction(next(value for end, value in CLIPPING if prime <= end), 100)
        upper, polynomial = extend(prime, clipping, Fraction(0), upper, polynomial)
        count += 1
    low, low_poly = upper, polynomial
    high, high_poly = upper, polynomial
    for prime in (7, 5):
        low, low_poly = extend(prime, Fraction(0), Fraction(0), low, low_poly)
        high, high_poly = extend(prime, Fraction(0), Fraction(prime - 1, prime), high, high_poly)
    positive = positive_multipliers(1, 3, high, high_poly, upper_weights(3))
    total = ceil(expectation(low[1], positive, Fraction(6, 5), 3))
    if count != 896 or total != 986888733243508975907980:
        raise ValueError("unexpected exact replay values")
    if 1000 * total >= 987 * SCALE:
        raise ValueError("reserve smaller than 13/1000")
    print(f"{count} large-prefix steps give total {total}/{SCALE} and reserve above 13/1000")
    verify_phase_constraints(Fraction(total, SCALE))


def verify_phase_constraints(base):
    root_cap = Fraction(4, 15)
    debit_45 = Fraction(8, 285)
    debit_63 = Fraction(4, 205)
    bound_without_15 = base + Fraction(6, 41) * root_cap - debit_45
    bound_without_21 = base + Fraction(4, 19) * root_cap - debit_45 - debit_63
    if bound_without_15 >= 1 - Fraction(2, 1000):
        raise ValueError("insufficient reserve when modulus 15 is absent")
    if bound_without_21 >= 1 - Fraction(4, 1000):
        raise ValueError("insufficient reserve when modulus 21 is absent")
    print("The phase constraints have reserve above 2/1000 and 4/1000")


if __name__ == '__main__':
    run()
