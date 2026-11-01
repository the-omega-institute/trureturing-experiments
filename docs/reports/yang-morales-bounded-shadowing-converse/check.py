#!/usr/bin/env python3
"""Exact checks for Yang–Morales, arXiv:2608.12862v1, Definitions 1 and 3.

These checks test finite instances and algebraic certificates of the proposed
argument. They do not decide its quantifiers over continuous tolerances or
over infinite pseudo-orbits. No floating point arithmetic is used.
"""

from fractions import Fraction as F
from itertools import product


def cadd(z, w):
    return z[0] + w[0], z[1] + w[1]


def cmul(z, w):
    return z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0]


def cscale(c, z):
    return c * z[0], c * z[1]


def norm_squared(z):
    return z[0] ** 2 + z[1] ** 2


def fidelity_isometry():
    """Lemma 5's ramp obstruction for four genuine complex isometries.

    On the radius-3 disk, delta(z)=1/(1+|z|^2) >= 1/10.
    A ramp from 0 to 3 in 31 steps has error 3/31 < 1/10.
    Epsilon=1 tracing is impossible: rotating both endpoint inequalities
    back would put the same point within 1 of both 0 and 3.
    """
    radius, steps = F(3), 31
    for lam in [(F(1), F(0)), (F(-1), F(0)),
                (F(0), F(1)), (F(0), F(-1))]:
        powers = [(F(1), F(0))]
        for _ in range(steps + 1):
            powers.append(cmul(lam, powers[-1]))
        orbit = [cscale(radius * F(n, steps), powers[n])
                 for n in range(steps + 1)]
        for n in range(steps):
            gx = cmul(lam, orbit[n])
            err = cadd(orbit[n + 1], cscale(F(-1), gx))
            delta_at_gx = F(1) / (1 + norm_squared(gx))
            assert norm_squared(err) == (radius / steps) ** 2
            assert norm_squared(err) < delta_at_gx ** 2
        assert norm_squared(orbit[0]) == 0
        assert norm_squared(orbit[-1]) == radius ** 2
        assert radius > 1 + 1
    print("FIDELITY: source Lemma 5, four unit-modulus rotations; "
          "strict ramp errors and incompatible endpoint balls verified")


def iterate(a, b, z, n):
    return a ** n * z[0], b ** n * z[1]


def norm_max(z):
    return max(abs(z[0]), abs(z[1]))


def distance(z, w):
    return norm_max((z[0] - w[0], z[1] - w[1]))


def bounded_recurrences(a, b):
    """All error words over {-delta/2,0,delta/2}, periods 1 through 3.

    The displayed cyclic sums are the exact infinite geometric-series
    values for these periodic errors. Check recurrence, bounds, and
    endpoint identities in the past for u and in the future for v.
    """
    delta = min(1 - a, b - 1) / 2  # epsilon(0)=1
    assert delta > 0
    alphabet = (-delta / 2, F(0), delta / 2)
    count = 0
    for period in range(1, 4):
        for errors in product(alphabet, repeat=period):
            def u(n):
                return sum((a ** j * errors[(n - 1 - j) % period]
                            for j in range(period)), F(0)) / (1 - a ** period)

            def v(n):
                return -sum((b ** (-j) * errors[(n + j - 1) % period]
                             for j in range(1, period + 1)), F(0)) / (1 - b ** (-period))

            for n in range(period):
                assert abs(errors[n]) < delta
                assert u(n + 1) == a * u(n) + errors[n]
                assert v(n + 1) == b * v(n) + errors[n]
                assert abs(u(n)) <= delta / (1 - a) <= F(1, 2)
                assert abs(v(n)) <= delta / (b - 1) <= F(1, 2)
                for length in (1, 2, 3, 7):
                    past = a ** length * u(n - length) + sum(
                        (a ** j * errors[(n - 1 - j) % period]
                         for j in range(length)), F(0))
                    future = b ** (-length) * v(n + length) - sum(
                        (b ** (-j) * errors[(n + j - 1) % period]
                         for j in range(1, length + 1)), F(0))
                    assert u(n) == past
                    assert v(n) == future
            count += 1
    # The two coordinate words are independent; their coordinate bounds
    # therefore also bound every joint two-coordinate periodic realization.
    return count


def splice_obstruction(a, b, scale, delta_power):
    """Variable delta(z)=scale/(1+||z||^h); source evaluation points.

    Negative-time exact orbit forces the stable coordinate of any tracer
    to be 1+eta; positive-time exact orbit forces its unstable coordinate
    to be 1 (since epsilon<=1). The program checks the splice's only error
    at n=-1, then a finite exact failure certificate for this forced tracer.
    The tail-forcing argument itself is supplied in the mathematical draft.
    """
    k = 1
    while a * b ** k <= 1:
        k += 1
    p = F(1), F(1)
    eta = F(1, 2)
    while not eta < scale / (1 + (1 + eta) ** delta_power):
        eta /= 2
    q = 1 + eta, F(1)

    def delta(z):
        return scale / (1 + norm_max(z) ** delta_power)

    def epsilon(z):
        return F(1) / (1 + norm_max(z) ** k)

    def x(n):
        return iterate(a, b, q if n < 0 else p, n)

    for n in range(-16, 17):
        gx = iterate(a, b, x(n), 1)
        error = distance(gx, x(n + 1))
        assert error == (eta if n == -1 else 0)
        assert error < delta(gx)
        if n == -1:
            assert gx == q  # delta is evaluated here, not at p.
    base = a * b ** k
    assert base > 1
    n = 0
    while eta * base ** n < 1:
        n += 1
    tracer_n = iterate(a, b, q, n)
    error = distance(tracer_n, x(n))
    assert error == eta * a ** n
    assert epsilon(tracer_n) < b ** (-k * n)
    assert error >= epsilon(tracer_n)
    # Exact witnesses to the two divergent tail factors used to force z=q.
    assert a ** (-1) > 1 and b > 1
    return k, n


def main():
    fidelity_isometry()
    stable = sorted({F(p, q) for q in range(2, 7) for p in range(1, q)})
    unstable = sorted({1 + F(p, q) for q in range(1, 6)
                       for p in range(1, q + 1)} | {F(3), F(4)})
    parameters = list(product(stable, unstable))
    recurrence_count = 0
    splice_count = 0
    largest_k = 0
    largest_n = 0
    for a, b in parameters:
        assert 0 < a < 1 < b
        recurrence_count += bounded_recurrences(a, b)
        for scale, h in [(F(1), 1), (F(1, 8), 2), (F(1, 32), 3)]:
            k, n = splice_obstruction(a, b, scale, h)
            largest_k, largest_n = max(largest_k, k), max(largest_n, n)
            splice_count += 1
    print(f"FAMILY_GRID: {len(stable)} stable rates x {len(unstable)} unstable rates "
          f"= {len(parameters)} nondegenerate pairs")
    print(f"BOUNDED: {recurrence_count} periodic coordinate error words; "
          "both recurrence directions, four endpoint lengths, zero-tracer bound <= 1/2")
    print(f"FULL: {splice_count} exact splice certificates for three variable delta laws; "
          f"largest selected k={largest_k}, largest failure time={largest_n}")
    print("SCOPE: finite exact checks of the proposed family mechanism; "
          "no universal tolerance, infinite-orbit, or Lean proof is claimed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
