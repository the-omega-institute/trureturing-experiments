#!/usr/bin/env python3
"""Exact checks for Lecomte's geometric-transition endpoint candidate.

These certificates and finite checks support a preregistration outline; they
do not compute C or C1 and do not prove the infinite geometric classification.
"""

from fractions import Fraction as Q
from itertools import combinations, product
from collections import Counter
import sympy as S


def potential_squared(d2):
    assert d2 > 0
    return Q(1) / d2**6 - Q(2) / d2**3


def distance_squared(p, q):
    return sum(((a - b)**2 for a, b in zip(p, q)), Q(0))


def energy(points):
    assert len(set(points)) == len(points)
    return sum((potential_squared(distance_squared(p, q))
                for p, q in combinations(points, 2)), Q(0))


def compositions(n):
    if n == 0:
        yield ()
    for k in range(1, n + 1):
        for tail in compositions(n - k):
            yield (k,) + tail


def layer_energy(profile, gamma, slope):
    return sum((gamma - slope * (n - 1) for n in profile), Q(0))


def main():
    # Source fidelity: geometric endpoint proposition, axial chain and κ*.
    chain3 = [(Q(1, 2) + j, Q(0), Q(0)) for j in range(3)]
    assert energy(chain3) == -Q(8319, 4096)
    assert (3 * -Q(1, 2) + Q(8319, 4096)) / 3 == Q(725, 4096)
    # Source's chain-with-arm witness: check geometry and energy at n=2,...,16.
    for n in range(2, 17):
        m = 2 * n
        points = [(Q(j) - Q(1, 2), Q(0), Q(0)) for j in range(1, m + 1)]
        points += [(Q(n) - Q(1, 2), Q(i), Q(0)) for i in range(1, n)]
        populations = Counter(p[0].numerator // p[0].denominator for p in points)
        assert sorted(populations) == list(range(m))
        assert populations[n - 1] == n
        assert len(points) == 3 * n - 1
        assert max(distance_squared(a, b) for a, b in combinations(points, 2)) == (m - 1)**2
        assert energy(points) <= -(3 * n - 2)
    print('FIDELITY: E(chain3)=-8319/4096; kappa_star(-1/2)=725/4096; arm witnesses n=2..16 OK')

    # Octahedron coordinates below are scaled by 1/sqrt(2); work with squared
    # distances divided by 2, avoiding any numerical square roots.
    vertices = [tuple(s if i == j else 0 for i in range(3))
                for j in range(3) for s in (-1, 1)]
    ds = [distance_squared(a, b) / 2 for a, b in combinations(vertices, 2)]
    assert Counter(ds) == {Q(1): 12, Q(2): 3}
    octahedron_energy = sum(map(potential_squared, ds), Q(0))
    assert octahedron_energy == -Q(813, 64)
    c_lower = -octahedron_energy / 6
    assert c_lower == Q(271, 128)
    # Integral-test tail: sum_{n>=3} n^-6 <= integral_2^infty x^-6 dx.
    zeta6_upper = 1 + Q(1, 64) + Q(1, 5 * 2**5)
    assert zeta6_upper == Q(327, 320)
    # zeta(12)>=1, so q=1+2*zeta(6)-zeta(12)<=2*zeta6_upper.
    q_upper = 2 * zeta6_upper
    k_upper = 1 + 2 * q_upper
    gap_lower = c_lower - q_upper
    deficit_lower = gap_lower / (2 * (gap_lower + k_upper))
    assert (q_upper, k_upper, gap_lower, deficit_lower) == (
        Q(327, 160), Q(407, 80), Q(47, 640), Q(47, 6606))
    assert deficit_lower > Q(1, 141)
    print('ENDPOINT CERTIFICATES: octahedron=-813/64; C>=271/128; q<=327/160; K<=407/80')
    print('GAP: C-C1>=47/640; eventual endpoint deficit fraction>=47/6606>1/141')

    # A self-contained coarse stability constant for the outline: cubes of side
    # 1/16 have within-cube distances <1/8. The negative inter-cube row sum
    # is bounded by W, using 24*m^2+2 lattice sites in shell m and zeta(4)<=4/3.
    intra = Q(8**12 - 2 * 8**6)
    row = Q(26) + Q(784, 3) * 16**6
    assert intra > row > 0
    for n in range(1, 101):
        assert 24 * (n + 1)**2 + 2 <= 98 * n**2
    print(f'STABILITY CERTIFICATE: intra={intra}; negative_row<={row}; Cstab={intra/2}')

    # Parameter identities for the entire affine zero-kernel family, including
    # gamma=-C1: these are symbolic exact polynomial identities.
    c, c1, k, gamma, n, d, u, delta = S.symbols('C C1 K gamma N D u delta', real=True)
    lb = gamma * d - (c + k) * (n - d)
    geometric_lb = -c1 * d - (c + k) * (n - d)
    assert S.expand(geometric_lb - lb - (-c1 - gamma) * d) == 0
    assert S.expand(lb + u * n - ((c + k + gamma) * d - (c + k - u) * n)) == 0
    threshold = (c + k - u) / (c + k + gamma)
    assert S.factor(1 - threshold - (c + gamma - delta) / (c + k + gamma)).subs(u, c-delta) == 0
    assert S.simplify(((c + gamma) / (2 * (c + k + gamma))).subs(gamma, -c1)
                      - (c-c1) / (2 * (c-c1+k))) == 0
    print('FAMILY: symbolic validity, strict exclusion, and endpoint deficit identities OK')

    # All ordered population profiles of total size 1..12, for nine parameter
    # choices and three exact-reference errors. These check strictness at equality.
    profiles_checked = 0
    equality_cases = 0
    for cc, ss in product((Q(271, 128), Q(3), Q(10)), (Q(1), Q(2), q_upper)):
        kk = k_upper
        for gg in (-ss, (-cc - ss) / 2, (-2 * cc - ss) / 3):
            beta = cc + gg
            assert beta > 0
            slope = cc + kk
            span_coefficient = slope + gg
            assert span_coefficient > 0
            for err in (Q(0), beta / 4, beta / 2):
                uu = cc - err
                for size in range(1, 13):
                    max_survivor = 0
                    for p in compositions(size):
                        span = len(p)
                        value = layer_energy(p, gg, slope)
                        assert value == gg * span - slope * (size-span)
                        assert value <= -ss * span - slope * (size-span)
                        excluded = value > -uu * size
                        assert excluded == (span_coefficient * span > (slope-uu) * size)
                        if value == -uu * size:
                            equality_cases += 1
                            assert not excluded
                        if not excluded:
                            max_survivor = max(max_survivor, span)
                        profiles_checked += 1
                    frac = beta / (2 * span_coefficient)
                    # floor threshold <= threshold; rho <= threshold whenever
                    # threshold>=0. The full bound permits zero surviving spans.
                    assert Q(size-max_survivor, size) >= frac
    # An explicit boundary sample checks the strict comparison independently
    # of whether the parameter grid happens to hit a rational equality.
    assert layer_energy((2,), Q(-1), Q(4)) == -Q(5, 2) * 2
    assert not (layer_energy((2,), Q(-1), Q(4)) > -Q(5, 2) * 2)
    equality_cases += 1
    print(f'PROFILE CHECK: {profiles_checked} exact cases; {equality_cases} equality cases remain nonexcluded')

    # Small actual oriented configurations: mandatory axial endpoints and
    # optional half-unit arms. Retain endpoints and one point per occupied layer.
    geometry_checked = 0
    for span in range(2, 6):
        endpoints = [(Q(1, 2), Q(0), Q(0)), (Q(span)-Q(1, 2), Q(0), Q(0))]
        middle = [(Q(i)+Q(1, 2), Q(0), Q(0)) for i in range(1, span-1)]
        retained = endpoints[:1] + middle + endpoints[1:]
        optional = [(Q(i)+Q(1, 2), Q(1, 2), Q(0)) for i in range(span)]
        for mask in range(1 << span):
            extras = [p for i, p in enumerate(optional) if (mask >> i) & 1]
            points = retained + extras
            diameter2 = max(distance_squared(a,b) for a,b in combinations(points,2))
            if diameter2 != (span-1)**2:
                continue
            assert energy(retained) >= -q_upper * span
            cross = sum((potential_squared(distance_squared(a,b)) for a in retained for b in extras), Q(0))
            assert cross >= -k_upper * len(extras)
            # This uses the certified coarse universal stability bound, not an
            # asserted numerical value of the unknown optimal constant C.
            assert energy(extras) >= -(intra/2) * len(extras)
            geometry_checked += 1
    print(f'GEOMETRY CHECK: {geometry_checked} oriented examples retain axial endpoints; singleton and cross-tail bounds OK')

    # Arbitrary nonzero-kernel two-pile algebra. For e=2t surplus particles,
    # populations are t+1 each; every other layer is a singleton.
    A, t, kap, D, G = S.symbols('A t kappa D G', real=True)
    bound = (D-2)*G + 2*A*(t+1) - kap*(t+1)**2
    expanded = G*D - 2*G + 2*A*t + 2*A - kap*t*t - 2*kap*t - kap
    assert S.expand(bound-expanded) == 0
    assert S.Poly(expanded, t).coeff_monomial(t*t) == -kap
    print('NONZERO KERNEL: two-pile quadratic coefficient is exactly -kappa(r), with linear pile costs')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
