#!/usr/bin/env python3
"""Exact checks for the Doric--Katic--Peric closing Question.

Only integers and fractions are used. Finite checks support, but do not prove,
the proposed infinite continuum, its covering dimension, or its entropy limit.
Run with Python 3; no external packages are required.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import gcd, isqrt


def require(test, message):
    if not test:
        raise AssertionError(message)


def ceilq(x):
    return -(-x.numerator // x.denominator)


def subsets(size):
    return range(1 << size)


def members(mask, size):
    return [i for i in range(size) if mask & (1 << i)]


def cover_number(masks, size):
    """Exact minimum of a cover, by dynamic programming on subsets."""
    full = (1 << size) - 1
    best = [size + 1] * (full + 1)
    best[0] = 0
    for old in subsets(size):
        for cover in masks:
            new = old | cover
            best[new] = min(best[new], best[old] + 1)
    return best[full]


def fidelity():
    """Source Section 1: periodic homeomorphisms have zero polynomial entropy.

    A non-isometric finite metric model has cycles of lengths 2 and 3. The
    exact period certificate implies d_N=d_6 for every N>=6. Therefore every
    separated cardinality stabilizes, reproducing this source result here.
    We also reproduce the source Section 2 cov/span/sep inequalities exactly.
    """
    size = 5
    f = [1, 0, 3, 4, 2]
    orbit = [[i] for i in range(size)]
    for i in range(size):
        for _ in range(18):
            orbit[i].append(f[orbit[i][-1]])
        require(orbit[i][6] == i, 'period-six certificate')
    d = lambda n, x, y: max(abs(orbit[x][j] - orbit[y][j]) for j in range(n))
    base = [[d(6, x, y) for y in range(size)] for x in range(size)]
    require(d(1, 2, 3) < d(6, 2, 3), 'fidelity model must be non-isometric')
    for n in range(6, 19):
        require([[d(n, x, y) for y in range(size)] for x in range(size)] == base,
                'Bowen distance stabilized after a full period')
    checks = 0
    last_sep = None
    for n in range(1, 19):
        for eps in [Q(1, 2), Q(1), Q(3, 2), Q(2), Q(5, 2)]:
            good = lambda mask, e: all(d(n, x, y) < e
                    for x, y in combinations(members(mask, size), 2))
            cov = cover_number([a for a in subsets(size) if good(a, eps)], size)
            cov2 = cover_number([a for a in subsets(size) if good(a, 2 * eps)], size)
            balls = [sum(1 << y for y in range(size) if d(n, x, y) < eps)
                     for x in range(size)]
            span = cover_number(balls, size)
            sep = max(len(members(a, size)) for a in subsets(size)
                      if all(d(n, x, y) >= eps
                             for x, y in combinations(members(a, size), 2)))
            require(cov2 <= span <= sep <= cov, 'source cov/span/sep chain')
            if n >= 6 and eps == Q(3, 2):
                require(sep == 4, 'stabilized separation cardinality')
                last_sep = sep
            checks += 1
    print('FIDELITY: source periodic zero-growth result; period=6; '
          f'sep(N,3/2)={last_sep} for N>=6 (exact period certificate).')
    print(f'FIDELITY: {checks} exact cov(2eps)<=span(eps)<=sep(eps)<=cov(eps) checks.')


def beta(n, k):
    return Q(1, 2 ** n) * (1 + Q(k, n))


def slope_checks():
    slopes = []
    for n in range(2, 129):
        for k in range(n):
            b = beta(n, k)
            require(Q(1, 2 ** n) <= b < Q(1, 2 ** (n - 1)), 'slope block')
            require(0 < b, 'positive slope')
            slopes.append(b)
    require(len(slopes) == len(set(slopes)), 'pairwise distinct spokes')
    for eps in [Q(1, 16), Q(1, 4), Q(1), Q(3)]:
        cutoff = 2
        while 6 * Q(1, 2 ** cutoff) >= eps / 8:
            cutoff += 1
        for n in range(cutoff, cutoff + 129):
            require(6 * Q(1, 2 ** n) < eps / 8, 'uniform tail height')
    print(f'GEOMETRY: {len(slopes)} distinct exact rational slopes; '
          'four uniform tail-height cutoffs checked.')


def modular_witness(n, m):
    g = gcd(n, m)
    q = m // g
    target = q // 2
    ell = (pow(n // g, -1, q) * target) % q
    j = n * ell
    require(1 <= ell < q, 'nonzero modular witness')
    require(j % n == 0 and (j % m) == g * target, 'orbit residues')
    require(0 <= j < n * m, 'strict observation-horizon bound')
    phase = Q(j % m, m)
    require(Q(1, 3) <= phase <= Q(1, 2), 'cosine nonpositive phase band')
    return j


def lower_checks():
    total = 0
    for n, m in combinations(range(2, 129), 2):
        modular_witness(n, m)
        total += 1
    # Actual Euclidean squared-distance checks for rational-cosine periods.
    # This table is exact, not a numerical evaluation of trigonometric functions.
    cosine = {Q(0): Q(1), Q(1, 6): Q(1, 2), Q(1, 4): Q(0),
              Q(1, 3): Q(-1, 2), Q(1, 2): Q(-1), Q(2, 3): Q(-1, 2),
              Q(3, 4): Q(0), Q(5, 6): Q(1, 2)}
    geometric = 0
    for c, r in [(Q(2), Q(1)), (Q(3), Q(1, 2)), (Q(7, 4), Q(5, 4))]:
        require(c > r > 0, 'nondegenerate family parameters')
        for n, m in combinations([2, 3, 4, 6], 2):
            j = modular_witness(n, m)
            xn = c + r * cosine[Q(j % n, n)]
            xm = c + r * cosine[Q(j % m, m)]
            yn = xn * beta(n, j % n)
            ym = xm * beta(m, j % m)
            require(xn - xm >= r, 'horizontal separation threshold r')
            require((xn - xm) ** 2 + (yn - ym) ** 2 >= r ** 2,
                    'Euclidean separation threshold r')
            geometric += 1
    for N in list(range(4, 1001)) + [4095, 4096, 4097, 65535, 65536, 65537]:
        M = isqrt(N)
        for n, m in combinations(range(2, M + 1), 2):
            j = modular_witness(n, m)
            require(j < n * m <= N, 'lower separated-set horizon')
    print(f'LOWER: {total} modular witnesses (2<=n<m<=128); '
          f'{geometric} exact Euclidean checks across three c>r>0 choices.')
    print('LOWER: floor(sqrt(N))-1 endpoint separation certificates for '
          '4<=N<=1000 and six square-boundary horizons.')


def upper_checks():
    horizons = list(range(2, 65)) + [127, 128, 255, 256, 257, 1023, 1024,
                                    1025, 4095, 4096, 4097, 65535, 65536,
                                    65537, 1000000]
    cells_checked = 0
    for eps in [Q(1, 16), Q(1, 4), Q(1), Q(3)]:
        for N in horizons:
            q = isqrt(N)
            ds, dt, dw = eps / 24, eps / 128, eps / (128 * N)
            # 2*pi<8 gives an exact rational overestimate for the analytic
            # Lipschitz bound. Analytic trigonometric facts are not inferred
            # from samples by this program.
            dx = 3 * ds + 8 * dt + 8 * N * dw
            dy = eps / 4
            require(dx == eps / 4 and dx ** 2 + dy ** 2 < eps ** 2,
                    'Bowen cell diameter envelope')
            # Keep each denominator <=q separately. All others lie in
            # [0,1/(q+1)]. No floating-point square root is used.
            tail_bins = ceilq(Q(1, q + 1) / dw) + 1
            bound = ceilq(128 * Q(q + 1) / eps) + 1
            require(N < (q + 1) ** 2 and tail_bins <= bound, 'reciprocal tail bound')
            for n in range(q + 1, q + 1025):
                index = Q(1, n) // dw
                require(0 <= index < tail_bins, 'tail frequency is covered')
            # Number of (s,theta) cells is fixed once epsilon is fixed.
            parameter_factor = (ceilq(1 / ds) + 1) * (ceilq(1 / dt) + 1)
            actual_bound = parameter_factor * (q + tail_bins)
            require(actual_bound <= parameter_factor * (q + bound), 'O_epsilon(sqrt(N)) count')
            cells_checked += 1
    print(f'UPPER: {cells_checked} rational cell-diameter and reciprocal-cover counts; '
          '1024 tail frequencies checked per case, N<=1000000.')


def main():
    fidelity()
    slope_checks()
    lower_checks()
    upper_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
