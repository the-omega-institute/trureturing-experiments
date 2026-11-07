#!/usr/bin/env python3
"""Exact finite-memory checks for Bertoin--Minchev, arXiv:2605.29030v1.

Only integer and Fraction arithmetic is used. Positive-vector comparisons
certify intervals for the spectral radius of the directly constructed transfer
matrix. An independent rational entropy value certifies a lower bound for B.
These finite checks are evidence for the proposed argument, not its proof.
"""
from fractions import Fraction as Q
from itertools import product
from math import lcm


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for k in range(total + 1):
            for rest in compositions(total - k, length - 1):
                yield (k,) + rest


def laws():
    ans = set()
    for length in range(1, 4):
        for denominator in range(1, 5):
            for v in compositions(denominator, length):
                if v[-1]:
                    ans.add(tuple(Q(x, denominator) for x in v))
    return sorted(ans, key=lambda t: (len(t), t))


def transfer(sigma, tau):
    states = list(product(range(len(sigma)), repeat=len(tau)))
    index = {s: i for i, s in enumerate(states)}
    matrix = [[Q(0) for _ in states] for _ in states]
    for i, s in enumerate(states):
        for t in range(len(sigma)):
            j = index[(t,) + s[:-1]]
            matrix[i][j] += sum(tau[k] * sigma[s[k]][t] for k in range(len(tau)))
        assert sum(matrix[i]) <= 1
    return matrix


def matvec(matrix, v):
    return [sum(a * x for a, x in zip(row, v)) for row in matrix]


def radius_interval(matrix, iterations=160):
    v = [Q(1)] * len(matrix)
    lower, upper = Q(0), Q(1)
    for _ in range(iterations):
        mv = matvec(matrix, v)
        assert min(v) > 0 and min(mv) > 0
        ratios = [a / b for a, b in zip(mv, v)]
        lower, upper = max(lower, min(ratios)), min(upper, max(ratios))
        v = mv
    assert lower <= upper
    return lower, upper


def entropy_power(sigma, p):
    """Return D, R with exp(D*F(p)) = R, exactly, for rational positive p."""
    q = [sum(p[i] * sigma[i][j] for i in range(len(p))) for j in range(len(p))]
    assert min(p) > 0 and min(q) > 0 and sum(p) == 1
    denominator = lcm(*(x.denominator for x in p))
    value = Q(1)
    a = [p[i] / q[i] for i in range(len(p))]
    for j in range(len(p)):
        assert sum(p[i] * sigma[i][j] * a[j] for i in range(len(p))) == p[j]
        value *= (q[j] / p[j]) ** int(p[j] * denominator)
    # The positive normalized LEFT eigenvector p has eigenvalue one.
    assert sum(p[i] * sum(sigma[i][j] * a[j] for j in range(len(p)))
               for i in range(len(p))) == 1
    return denominator, value


def entropy_certificate(sigma, upper):
    for d in range(2, 41):
        for k in range(1, d):
            p = (Q(k, d), Q(d - k, d))
            denominator, value = entropy_power(sigma, p)
            if upper ** denominator <= value:
                return p, denominator, value
    raise AssertionError("No exact ceiling certificate on the specified rational grid")


def above_benchmark(x, sigma):
    # For a positive off-diagonal 2x2 matrix: r=(trace+sqrt(discriminant))/2.
    a, b = sigma[0]
    c, d = sigma[1]
    u = 2 * x - a - d
    return u > 0 and u * u > (a - d) ** 2 + 4 * b * c


def display_interval(interval):
    scale = 10 ** 9
    lo, hi = interval
    return (f"[{lo.numerator * scale // lo.denominator}/{scale}, "
            f"{(hi.numerator * scale + hi.denominator - 1) // hi.denominator}/{scale}]")


def boundary_check(sigma):
    count = 0
    for vertex in range(2):
        if sigma[vertex][vertex] == 0:
            # This boundary entropy is -infinity, not Real.log(0)=0.
            continue
        epsilon = Q(1, 100)
        p = [epsilon, epsilon]
        p[vertex] = 1 - epsilon
        denominator, value = entropy_power(sigma, p)
        assert value > sigma[vertex][vertex] ** denominator
        count += 1
    return count


def telescoping_checks(tau):
    b = (Q(-2), Q(3))
    count = 0
    for N in range(4):
        width = max(N + 1, len(tau))
        w = lambda j: sum(tau[j:], Q(0))
        for s in product(range(2), repeat=width):
            z = sum(tau[j] * b[s[j]] for j in range(len(tau)))
            delta = sum((tau[j] * (b[s[j]] - b[s[N]])
                         for j in range(N, len(tau))), Q(0))
            exponent = sum(w(j) * b[s[j]] for j in range(N + 1))
            assert delta <= (max(b) - min(b)) * w(N)
            for t in range(2):
                prefixed = (t,) + s
                shifted = sum(w(j) * b[prefixed[j]] for j in range(N + 1))
                assert shifted - exponent == b[t] - z + delta
                count += 1
    return count


def main():
    matrices = {
        "source_figures_1_2": ((Q(18, 25), Q(2, 25)), (Q(9, 50), Q(29, 50))),
        "positive_asymmetric": ((Q(3, 5), Q(1, 5)), (Q(1, 10), Q(3, 5))),
        "irreducible_with_zero": ((Q(0), Q(1, 2)), (Q(1, 3), Q(1, 4))),
    }
    # Source's displayed matrix has benchmark r=(65+sqrt(193))/100.
    source = matrices["source_figures_1_2"]
    trace = source[0][0] + source[1][1]
    determinant = source[0][0] * source[1][1] - source[0][1] * source[1][0]
    assert trace == Q(13, 10) and determinant == Q(1008, 2500)
    assert trace * trace - 4 * determinant == Q(193, 2500)
    assert transfer(source, (Q(1),)) == [list(row) for row in source]
    print("FIDELITY source matrix: benchmark r=(65+sqrt(193))/100; tau=delta_0 gives exactly sigma")
    samples = laws()
    strict_count = 0
    boundary_count = 0
    representatives = [(Q(1),), (Q(1, 2), Q(1, 2)), (Q(1, 3),) * 3]
    for name, sigma in matrices.items():
        assert sigma[0][1] > 0 and sigma[1][0] > 0
        assert any(sigma[i][i] > 0 for i in range(2))
        assert sum(sigma[0]) != sum(sigma[1])
        assert all(x >= 0 for row in sigma for x in row)
        boundary_count += boundary_check(sigma)
        for tau in samples:
            matrix = transfer(sigma, tau)
            interval = radius_interval(matrix)
            p, denominator, value = entropy_certificate(sigma, interval[1])
            assert interval[1] ** denominator <= value
            if sum(x > 0 for x in tau) >= 2:
                assert above_benchmark(interval[0], sigma)
                strict_count += 1
            if tau in representatives:
                print(f"EXACT {name}: tau={tuple(map(str,tau))}; dimension={len(matrix)}; "
                      f"radius in {display_interval(interval)}; "
                      f"exp({denominator}*F({tuple(map(str,p))}))={value}; upper^{denominator}<=certificate")
        print(f"CEILING {name}: all {len(samples)} laws of support within {{0,1,2}}, reduced denominator<=4 passed")
    telescopes = sum(telescoping_checks(tau) for tau in samples)
    print(f"FIDELITY Theorem 4.1(ii): strict enhancement certified in all {strict_count} non-Dirac cases")
    print(f"BOUNDARY {boundary_count} finite-entropy vertices improve under epsilon=1/100; zero-diagonal vertex has entropy -infinity")
    print(f"TELESCOPING {telescopes} exact prefix identities and tail-error inequalities passed, including N=0")
    print(f"SCOPE {3 * len(samples)} direct finite-memory spectral intervals; no universal conclusion follows from enumeration")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
