#!/usr/bin/env python3
"""Exact finite checks for Yang's Question 5 and its proposed counting argument.

Only integers and Fraction are used. These checks do not establish convergence.
Source: Kai Yang, Analysis of Convergence for Some Positive Oscillating Series
and Integrals, Memorial MSc thesis, Chapter 4 Lemma 4.1 and Chapter 6 Question 5.
"""
from fractions import Fraction as F
from math import gcd, isqrt


def circle_dist(x):
    r = x - x.numerator // x.denominator
    return min(r, 1 - r)


def fidelity_dirichlet():
    # Lemma 4.1: some 1 <= q <= N has ||q theta|| <= 1/(N+1).
    # theta=sqrt(2). Squaring rational endpoints decides every comparison exactly.
    witnesses = []
    for N in range(1, 513):
        t = F(1, N + 1)
        for q in range(1, N + 1):
            lo = isqrt(2 * q * q)
            m = lo if 8 * q * q < (2 * lo + 1) ** 2 else lo + 1
            if (m - t) ** 2 <= 2 * q * q <= (m + t) ** 2:
                witnesses.append(q)
                break
        else:
            raise AssertionError((N, "Dirichlet witness missing"))
    print("FIDELITY Lemma 4.1, theta=sqrt(2): N=1..512, 512 exact witnesses")


def kernel_deficit():
    # Proposed envelope: |(a+y)/(a+1)| <= 1-kappa(1-y), -1<=y<=1.
    # Includes 0<a<1, where omitting absolute values would be a defect.
    count = 0
    for a in [F(i, j) for i in range(1, 13) for j in range(1, 9)]:
        kappa = min(a, F(1)) / (a + 1)
        for j in range(-128, 129):
            y = F(j, 128)
            lhs = abs((a + y) / (a + 1))
            rhs = 1 - kappa * (1 - y)
            assert 0 <= lhs <= rhs <= 1, (a, y)
            for n in (1, 2, 3, 8, 16):
                assert lhs ** n <= rhs ** n
            count += 1
    print(f"KERNEL deficit and integer powers: {count} exact (a,y) pairs")


def perturbed_orbit_counts():
    # Each length-q orbit block differs from a translated q-grid by <=1/Q.
    # Count a closed arc about beta=1/4, including the wrap at 0.
    # The bound 2*q*r+4 pays both grid endpoints and the perturbation.
    count = 0
    for q in range(1, 21):
        for p in range(1, q + 1):
            if gcd(p, q) != 1:
                continue
            for Q in sorted({q, 2 * q, q * q}):
                for z in (F(0), F(1, 4), F(1, 3 * q)):
                    for sign in (-1, 0, 1):
                        eps = F(sign, q * Q)
                        distances = [circle_dist(z + i * (F(p, q) + eps) - F(1, 4))
                                     for i in range(q)]
                        for radius in (F(0), F(1, 64), F(1, 16), F(1, 4), F(1, 2)):
                            hit = sum(d <= radius for d in distances)
                            assert hit <= 2 * q * radius + 4, (p, q, Q, z, sign, radius)
                            # Every prefix is a subset of this same block.
                            for length in (0, q // 2, q):
                                assert sum(d <= radius for d in distances[:length]) <= 2*q*radius+4
                            count += 1
    print(f"ORBIT closed-arc counts, all coprime 1<=p<=q<=20: {count} configurations")


def quadratic_phase_block_mass():
    # alpha=sqrt(2)-1, theta=2*pi*alpha, beta=1/4.
    # Rational enclosure of sqrt(2), certified by integer squares.
    B = 10 ** 30
    root = isqrt(2 * B * B)
    assert root * root < 2 * B * B < (root + 1) ** 2
    A = root - B
    D = 4 * B
    precision = 10 ** 9
    for k in range(1, 7):
        M = 4 ** k
        sqrtM = 2 ** k
        upper_units = 0
        for n in range(M, 2 * M):
            lo, hi = 4 * n * A - B, 4 * n * (A + 1) - B
            # A certified lower bound on ||n*alpha-beta||.
            candidates = range(lo // D, hi // D + 2)
            d = min(0 if lo <= m * D <= hi else min(abs(lo-m*D), abs(hi-m*D))
                    for m in candidates)
            # a=1/4, kappa=1/5, p=1. The Gaussian envelope is <=
            # 1/(1+8*kappa*M*d^2). Round each rational UPPER bound upwards.
            num = 5 * D * D * precision
            den = 5 * D * D + 8 * M * d * d
            upper_units += (num + den - 1) // den
        assert upper_units <= 10 * sqrtM * precision, (M, upper_units)
        print(f"BLOCK M={M}: mass majorant <= {upper_units}/{precision} <= 10*sqrt(M)={10*sqrtM}")


def main():
    fidelity_dirichlet()
    kernel_deficit()
    perturbed_orbit_counts()
    quadratic_phase_block_mass()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
