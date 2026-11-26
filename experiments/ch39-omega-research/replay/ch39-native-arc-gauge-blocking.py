"""Finite character blocking at fixed native gauge trace stiffness."""
from fractions import Fraction as Q

kappa = Q(1, 480)


def irrep(group, n):
    if group == 'U1':
        return Q(n*n), 2
    if group == 'SU2':
        return Q(n*(n+2), 3), (n+1)**2
    # Symmetric power (n,0) and its conjugate (0,n).
    dim = (n+1)*(n+2)//2
    return Q(n*(n+3), 4), 2*dim**2


def stiffness(sectors, cells):
    p0 = 1 - sum(p for _, _, p in sectors)
    assert p0 > 0 and all(p >= 0 for _, _, p in sectors)
    terms = [p**cells / m**(cells-1) for c, m, p in sectors]
    return sum(c*v for (c, _, _), v in zip(sectors, terms, strict=True))/(p0**cells + sum(terms))


def groups():
    for group, k, unit, sharp in [
        ('U1', 120*kappa, Q(1, 60), Q(1, 1140)),
        ('SU2', 8*kappa, Q(1, 4), Q(1, 55700)),
        ('SU3', 12*kappa, Q(1, 6), Q(1, 164274)),
    ]:
        c1, m1 = irrep(group, 1)
        c32, m32 = irrep(group, 32)
        low = [(c32, m32, k/c32)]
        high = [(c1, m1, k)]
        assert c1 == 1 and k*unit == Q(1, 240)
        assert stiffness(low, 1) == stiffness(high, 1) == k
        assert unit*stiffness(low, 2) < Q(1, 10**6) < Q(5, 10**6) < sharp
        assert unit*stiffness(high, 2) == sharp
        for cells in (2, 3, 8):
            opt = k**cells / (m1**(cells-1)*(1-k)**cells + k**cells)
            assert stiffness(high, cells) == opt
        print(group, 'two-cell upper endpoint', sharp, 'degree-32 witness', unit*stiffness(low, 2))


def finite_phases():
    # w_j = exp(i phi_j), with phases in eighths of a turn on four group elements.
    # Accumulate the convolution directly in the basis 1,z,z^2,z^3 of Z[z]/(z^4+1).
    def blocked(phases, output):
        coeff = [0]*4
        for j in range(4):
            n = (phases[j] + phases[(output-j) % 4]) % 8
            coeff[n % 4] += (1 if n < 4 else -1)
        return coeff
    for phases, identity, opposite in [
        ([0, 1, 2, 1], [0, 0, 2, 0], [0, 0, 4, 0]),
        ([0, 2, 6, 2], [-2, 0, 0, 0], [-2, 0, -2, 0]),
    ]:
        assert blocked(phases, 0) == identity
        assert blocked(phases, 2) == opposite
    print('Four holonomies, opposite-to-identity amplitude ratios 2 and 1+i.')


if __name__ == '__main__':
    groups()
    finite_phases()
