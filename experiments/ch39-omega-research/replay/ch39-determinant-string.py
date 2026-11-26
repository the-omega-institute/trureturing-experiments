"""Finite certificates for the determinant ring and the three mass pencils."""
from fractions import Fraction
from itertools import combinations
from math import prod
from random import Random

from flint import arb, fmpq, fmpq_mat, fmpq_poly


def check_ring():
    # Turns, with D = d - i A_t on the target. A flat source holonomy is allowed.
    weights = [Fraction(s, 37) for s in (2, -3, 7, 5, -11)]
    pairs = list(combinations(range(5), 2))
    lifted = [weights[i] + weights[j] for i, j in pairs]
    assert sum(lifted) == 0
    for k in (-13, -1, 0, 1, 13):
        winding = [k] + [0]*9
        target = [w - Fraction(k, 10) + s for w, s in zip(winding, lifted)]
        assert sum(target) == 0
        residual = [w - t + s for w, t, s in zip(winding, target, lifted)]
        assert sum(h*h for h in residual) == Fraction(k*k, 10)
        assert all((t-s) % 1 == Fraction(-k, 10) % 1 for t, s in zip(target, lifted))
    assert 20*(1-arb.cos_pi_fmpq(fmpq(1, 10))) < 1
    for n in (3, 4, 5, 6, 7, 9, 19):
        # Read the cube of the uniformly advancing determinant using shortest arcs.
        turns = Fraction(3, n) % 1
        if turns == Fraction(1, 2):
            assert n == 6
        else:
            observed = n*(turns if turns < Fraction(1, 2) else turns-1)
            assert observed == {3: 0, 4: -1, 5: -2, 7: 3, 9: 3, 19: 3}[n]
    print("Unit winding on 3, 6, 9 uniform edges: cubic readout 0, ambiguous, 3.")


def check_family_pencils():
    rng = Random(390177)
    n = 10
    x = [2, -3, 5]
    y = [7, 11, -13]
    ratio = fmpq(prod(x), prod(y))
    X, Y = [fmpq_mat([[rng.randrange(-4, 5) for _ in range(n)]
                       for _ in range(n)]) for _ in range(2)]
    assert X*Y != Y*X
    mass = fmpq_mat(3*n, 3*n)
    for a in range(3):
        for i in range(n):
            for j in range(n):
                mass[((-a-1) % 3)*n+i, a*n+j] = x[a]*X[i, j]
                mass[(-a % 3)*n+i, a*n+j] = y[a]*Y[i, j]
    # Interpolate the degree <=10 pencil, then take its norm over lambda^3 = ratio.
    points = range(n+1)
    vandermonde = fmpq_mat([[t**d for d in points] for t in points])
    coeff = vandermonde.solve(fmpq_mat([[(Y+t*X).det()] for t in points]))
    pencil = fmpq_poly([coeff[d, 0] for d in points])
    norm = fmpq_poly([-ratio, 0, 0, 1]).resultant(pencil)
    assert mass.det() == prod(y)**n*norm != 0
    # A separate vacuum example has rank P=3 and a genuinely mixed orthogonal U.
    P = fmpq_mat([[int(i == j and i < 3) for j in range(n)] for i in range(n)])
    I = fmpq_mat([[int(i == j) for j in range(n)] for i in range(n)])
    U = I - fmpq_mat([[fmpq(2, n)]*n for _ in range(n)])
    assert U.transpose()*U == I and U.det() == -1
    assert all((U*(I-P+t*P)).det() == U.det()*t**3 for t in points)
    print("Dense 30 by 30 determinant agrees exactly with the three-pencil resultant.")


def check_disk():
    # Write c=lambda*a at a split core. Each branch has |c|=2/3.
    # f=z+c(1-|z|^2) forces z=-r*c/|c| and has nonsingular real derivative.
    s, r = fmpq(2, 3), fmpq(1, 2)
    assert r == s*(1-r*r)
    assert 1 + 2*s*r == fmpq(5, 3)
    # The three core points solve z^3=-(a^3*ratio)*(1-r^2)^3.
    assert s**3*(1-r*r)**3 == r**3
    print("Three separated simple zeros at radius 1/2, each of local degree +1.")


if __name__ == "__main__":
    check_ring()
    check_family_pencils()
    check_disk()
