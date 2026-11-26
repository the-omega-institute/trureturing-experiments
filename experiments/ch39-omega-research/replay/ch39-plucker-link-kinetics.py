"""Finite-band phase transport on all pairs of weak planes and an exact spectator marginal."""
from fractions import Fraction as F
from itertools import product
from math import prod


def positive_path_minors(size, x):
    # Sylvester minors of xI minus the half-weight path adjacency.
    a, b = F(1), x
    for _ in range(1, size):
        a, b = b, x*b-a/4
        if b <= 0:
            return False
    return x > 0


def verify():
    # Two-plane stabilizers for orthogonal planes fix g while demanding opposite link phases.
    h = (1, 1, -1, -1, 1)
    assert prod(h) == 1 and h[3]*h[4] == -1 and h[0]*h[1] == 1
    c = list(range(1, 9))+list(range(8, 0, -1))
    scale = sum(a*a for a in c)
    first = F(sum(a*b for a, b in zip(c, c[1:])), scale)
    assert first == F(50, 51)
    # At bandwidth 36 all normalized positive kernels lose more than 1% on a triangle.
    # Bandwidth 37 can lose less than 1%.
    assert positive_path_minors(37, F('.99665')) and F('.99665')**3 < F('.99')
    assert not positive_path_minors(38, F('.99666')) and F('.99666')**3 > F('.99')

    # Independent finite Fourier check: a Z4 triangle with nonconstant relative plaquette weight.
    # Link kernel 1+Re(b*conj(eta))/2 sampled at the four roots; b is real here.
    cosine = (1, 0, -1, 0)
    overlaps = (F(1, 2), F(-1, 4), F(3, 4))
    relative = (F(3, 2), F(1), F(1, 2), F(1))
    total = mixed = F(0)
    for etas in product(range(4), repeat=3):
        conditional = prod(1+b*cosine[k]/2 for b, k in zip(overlaps, etas))
        for links in product(range(4), repeat=3):
            h = relative[(sum(links)-sum(etas)) % 4]
            total += conditional*h
            mixed += conditional*h*cosine[sum(links) % 4]
    assert total == 4**6
    relative_first = sum(a*b for a, b in zip(relative, cosine))/sum(relative)
    assert mixed/total == relative_first*prod(overlaps)/4**3
    print('Finite triangle identity:', mixed/total, '15-mode rational first moment:', first,
          'smallest bandwidth for 99% triangle retention:', 37)


if __name__ == '__main__':
    verify()
