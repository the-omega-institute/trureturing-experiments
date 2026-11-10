"""Exact fidelity and counterexample checks for Çevik 10.1–10.2."""
from fractions import Fraction as Q
from itertools import product
from math import prod


def check_pyramid():
    """Example 6.4: two image majorants with coordinatewise meet y."""
    for y in product(range(-2, 3), repeat=4):
        delta = y[0] + y[1] - y[2] - y[3]
        indices = (2, 3) if delta >= 0 else (0, 1)
        majorants = [
            tuple(yk + (abs(delta) if k == j else 0)
                  for k, yk in enumerate(y))
            for j in indices
        ]
        for z in majorants:
            a, b, c = Q(z[0] - z[1], 2), Q(z[2] - z[3], 2), Q(z[0] + z[1], 2)
            assert (c + a, c - a, c + b, c - b) == z
            assert all(zk >= yk for zk, yk in zip(z, y))
        assert tuple(map(min, zip(*majorants))) == y
    print("Source Example 6.4: 625 exact coordinatewise meet checks")


def spike(n, t):
    return max(1 - n * t, 0) + Q(1, n)


def check_spikes():
    for n in range(1, 129):
        # Both functions are affine between successive breakpoints.
        for t in {Q(0), Q(1, n + 1), Q(1, n), Q(1)}:
            assert Q(1, n + 1) <= spike(n + 1, t) <= spike(n, t)
        assert spike(n, Q(0)) == 1 + Q(1, n)
        assert spike(n, Q(1, n)) == Q(1, n)
    print("Spikes n=1..129: strict positivity, monotonicity, retained value at 0")


def check_quadratic_majorants():
    # The spike is n-Lipschitz. For epsilon>0, epsilon+M*s^2 >= n*|s|
    # for all s when 4*M*epsilon >= n^2, by completing the square.
    for n, j, den in product(range(1, 9), range(9), (2, 5, 17)):
        x, epsilon = Q(j, 8), Q(1, den)
        M = n * n / (4 * epsilon)
        for lo, hi, slope, offset in (
            (Q(0), Q(1, n), -n, 1 + Q(1, n)),
            (Q(1, n), Q(1), 0, Q(1, n)),
        ):
            def gap(t):
                return spike(n, x) + epsilon + M * (t - x)**2 - (slope*t + offset)

            vertex = x + slope / (2 * M)
            assert min(gap(lo), gap(hi), gap(min(hi, max(lo, vertex)))) >= 0
    print("Quadratic order-density witnesses: 216 exact continuous-interval checks")


def multiply_by_root(coefficients, root):
    out = [Q(0)] * (len(coefficients) + 1)
    for j, value in enumerate(coefficients):
        out[j] -= root * value
        out[j + 1] += value
    return out


def check_positive_tests():
    for m in range(2, 9):
        nodes = [Q(j, m + 1) for j in range(1, m + 2)]
        bounds = []
        for j, x in enumerate(nodes):
            others = [z for k, z in enumerate(nodes) if j != k]
            denominator = prod(x - z for z in others)
            polynomial = [Q(1)]
            for z in others:
                polynomial = multiply_by_root(polynomial, z)
            values = [
                sum(c * t**d for d, c in enumerate(polynomial)) / denominator
                for t in nodes
            ]
            assert values == [int(k == j) for k in range(len(nodes))]
            # For t in [0,1], |t-z| <= max(z, 1-z). Thus the sum bounds
            # the norm via exact interpolation without sampling the interval.
            bounds.append(prod(max(z, 1-z) for z in others) / abs(denominator))
        C = sum(bounds)
        assert C >= 1
        for n in (m + 1, 2*(m + 1), 8*(m + 1)):
            assert all(spike(n, t) == Q(1, n) for t in nodes)
        print(f"Degree {m}: verified interpolation bound ||e|| <= ({C})/n "
              f"for 0 <= e <= b_n, n >= {m + 1}")


if __name__ == "__main__":
    check_pyramid()
    check_spikes()
    check_quadratic_majorants()
    check_positive_tests()
    print("ALL_CHECKS_OK")
