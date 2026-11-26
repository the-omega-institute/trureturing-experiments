"""Rank-four Gaussian Wilson coefficient and uniform finite-strength error budgets."""
from fractions import Fraction as F


def covariance():
    # Relative to the Hilbert--Schmidt metric the precision eigenvalues of Q
    # are 1 on su(4), 1/2 on mixed directions and 1/5 on the trace direction.
    spectrum = [(15, F(1)), (8, F(1, 2)), (1, F(1, 5))]
    variances = [1/(2*precision) for count, precision in spectrum for _ in range(count)]
    second = sum(variances)
    fourth = second**2+2*sum(v*v for v in variances)
    assert second == 18 and fourth <= 972
    assert second/10 == F(9, 5)
    return second


def error_bounds():
    degree = 14
    n = 8
    # At larger n the left side grows fourfold per step, the budget is affine.
    tail_budget = 23*degree + 171 + 130*(n+1) + 3*n
    assert F(4**n, 20) > tail_budget
    assert 3*tail_budget > 133
    assert F(18, 4**n) < F(1, 2)
    assert 16*14 + 10 + F(2**10, 12*4**n) < 240
    assert 2*240*F(2**n, 2**(4*n)) < F(1, 2)
    assert 1 + 1024 + 260 + F(4*240*2, 4**n) < 4096
    # A sharper budget for a triangle at the published finite example.
    strength = 2**115
    radius = F(1, 2**46)
    scaled2 = strength*radius**2
    steps = 3
    tail = F(1, 2**30)
    assert scaled2/20 > 23*degree + 171 + 13*115 + 2 + 30
    error = (tail + strength*(steps*radius)**3
             + 260*steps**2/scaled2 + 4*240*steps*radius*2*steps**2)
    assert error < F(1, 2500)
    print(f"At kappa=2^115 the triangle coefficient error is below 1/2500 (budget {float(error):.7f}).")


if __name__ == "__main__":
    covariance()
    error_bounds()
