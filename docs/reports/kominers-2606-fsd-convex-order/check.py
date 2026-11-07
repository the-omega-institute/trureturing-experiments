#!/usr/bin/env python3
"""Exact checks for Kominers arXiv:2606.03482v2, Section 5's open extension.

Standard library only. All evidence uses integers and fractions.Fraction.
Finite checks support a preregistration; they are not a Lean settlement.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from math import comb


def span(rows):
    words = {0}
    for row in rows:
        words |= {x ^ row for x in words}
    return words


def dual(words, n):
    return {x for x in range(1 << n)
            if all((x & c).bit_count() % 2 == 0 for c in words)}


def spectrum(words, n):
    counts = Counter(x.bit_count() for x in words)
    return tuple(counts[w] for w in range(n + 1))


def convolution(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def polynomial(coeffs, z):
    return sum((a * z ** j for j, a in enumerate(coeffs)), Q(0))


def mean_variance(weight_counts, scale=1):
    size = sum(weight_counts)
    mean = sum(Q(w, scale) * a for w, a in enumerate(weight_counts)) / size
    var = sum((Q(w, scale) - mean) ** 2 * a
              for w, a in enumerate(weight_counts)) / size
    return mean, var


def convex_order(weight_counts, k):
    """Exact criterion: equal means and all hinge inequalities at joint knots.

    Compare K = weight/2 with B ~ Bin(k,1/2). Both laws have size 2^k.
    Every piecewise linear hinge gap is linear between consecutive half-integers.
    Thus checking these knots checks every hinge, hence every convex test on
    the finite joint support (linear interpolation supplies the criterion).
    """
    assert sum(weight_counts) == 1 << k
    benchmark = [comb(k, j) for j in range(k + 1)]
    if sum(w * a for w, a in enumerate(weight_counts)) != k * (1 << k):
        return False
    return all(sum(a * max(t - w, 0) for w, a in enumerate(weight_counts))
               <= sum(a * max(t - 2 * j, 0) for j, a in enumerate(benchmark))
               for t in range(2 * k + 1))


def fidelity_h8():
    # RM(1,3): evaluations of all affine binary functions on F_2^3.
    points = list(product(range(2), repeat=3))
    rows = [(1 << 8) - 1]
    rows += [sum(x[j] << i for i, x in enumerate(points)) for j in range(3)]
    words = span(rows)
    expected = (1, 0, 0, 0, 14, 0, 0, 0, 1)
    assert len(words) == 16 and words == dual(words, 8)
    assert spectrum(words, 8) == expected
    reference = (1, 0, 4, 0, 6, 0, 4, 0, 1)
    gap = tuple(a - b for a, b in zip(reference, expected))
    assert gap == (0, 0, 4, 0, -8, 0, 4, 0, 0)
    assert convex_order(expected, 4)
    print('FIDELITY H8: self-dual; W(1,z)=1+14*z^4+z^8; '
          'P(z)=4*z^2*(1-z^2)^2; all hinge tests pass')


def family():
    for m in range(3, 65):
        rep_counts = tuple(1 if w in (0, m) else 0 for w in range(m + 1))
        even_counts = tuple(comb(m, w) if w % 2 == 0 else 0
                            for w in range(m + 1))
        assert sum(even_counts) == 1 << (m - 1)
        counts = convolution(rep_counts, even_counts)
        assert sum(counts) == 1 << m
        assert counts == convolution(even_counts, rep_counts)
        mu_r, var_r = mean_variance(rep_counts)
        mu_e, var_e = mean_variance(even_counts)
        mu_k, var_k = mean_variance(counts, scale=2)
        mu_b, var_b = mean_variance(tuple(comb(m, j) for j in range(m + 1)))
        assert mu_r == mu_e == mu_k == mu_b == Q(m, 2)
        assert var_r == Q(m * m, 4) and var_e == Q(m, 4)
        assert var_k == Q(m * (m + 1), 16) and var_b == Q(m, 4)
        assert var_k - var_b == Q(m * (m - 3), 16)
        assert (var_k > var_b) == (m >= 4)
        # A parity condition leaves 2^(m-r-1) completions on r fixed bits.
        for r in (1, 2):
            for fixed in product(range(2), repeat=r):
                remaining = m - r
                count = sum(comb(remaining, w) for w in range(remaining + 1)
                            if (w + sum(fixed)) % 2 == 0)
                assert count == 1 << (m - r - 1)
        if m >= 4:
            assert not convex_order(counts, m)
    print('FAMILY m=3..64: exact coefficient laws, parity fibers, means, variances; '
          'variance excess=m*(m-3)/16; strict for every checked m>=4')

    # Independently enumerate the ambient orthogonal complement, rather than
    # obtain it from the proposed block-swap identity, for small m.
    for m in range(3, 9):
        mask = (1 << m) - 1
        rep = {0, mask}
        even = {x for x in range(1 << m) if x.bit_count() % 2 == 0}
        words = {r | (e << m) for r in rep for e in even}
        expected_dual = {e | (r << m) for e in even for r in rep}
        assert len(words) == 1 << m
        assert dual(words, 2 * m) == expected_dual
        swapped = {((x & mask) << m) | (x >> m) for x in expected_dual}
        assert swapped == words
        assert spectrum(words, 2 * m) == spectrum(expected_dual, 2 * m)
        assert min(x.bit_count() for x in words if x) == 2
        assert all(any(x & (1 << j) for x in words) for j in range(2 * m))
        assert all(x.bit_count() % 2 == 0 for x in words) == (m % 2 == 0)
        print(f'ENUMERATED m={m}: |C|={len(words)}; full support; d=2; '
              f'isodual; even={m % 2 == 0}')


def witness():
    rows = [int(s, 2) for s in ('11110000', '00001100', '00001010', '00001001')]
    words = span(rows)
    counts = spectrum(words, 8)
    assert len(words) == 16 and counts == (1, 0, 6, 0, 2, 0, 6, 0, 1)
    assert counts == spectrum(dual(words, 8), 8)
    assert words != dual(words, 8)
    assert all(x.bit_count() % 2 == 0 for x in words)
    mu, var = mean_variance(counts, scale=2)
    assert mu == 2 and var == Q(5, 4)
    assert mean_variance((1, 4, 6, 4, 1)) == (Q(2), Q(1))
    reference = (1, 0, 4, 0, 6, 0, 4, 0, 1)
    difference = tuple(a - b for a, b in zip(counts, reference))
    assert difference == (0, 0, 2, 0, -4, 0, 2, 0, 0)
    z = Q(1, 2)
    assert polynomial(counts, z) == Q(697, 256)
    assert (1 + z * z) ** 4 == Q(625, 256)
    assert polynomial(difference, z) == Q(9, 32)
    print('WITNESS words=' + ','.join(f'{x:08b}' for x in sorted(words)))
    print('WITNESS half-weight multiplicities=(1,6,2,6,1)/16; '
          'mean=2; E[(K-2)^2]=5/4 > E[(B-2)^2]=1')
    print('WITNESS W-(1+z^2)^4=2*z^2*(1-z^2)^2; '
          'at z=1/2: 697/256-625/256=9/32')


def strengthenings():
    # Exhaust all systematic [2k,k] binary generator presentations [I|A], k<=4.
    # Any such code can be made systematic by a coordinate permutation, which
    # preserves all the properties tested. Counts are presentations, not classes.
    for k in range(1, 5):
        mask = (1 << k) - 1
        count_fsd_d3 = count_doubly_even_fsd = 0
        failures = 0
        for encoded in range(1 << (k * k)):
            rows_a = [(encoded >> (i * k)) & mask for i in range(k)]
            rows_at = [sum(((rows_a[j] >> i) & 1) << j for j in range(k))
                       for i in range(k)]
            words = span([(1 << i) | (rows_a[i] << k) for i in range(k)])
            dual_words = span([rows_at[i] | (1 << (k + i)) for i in range(k)])
            counts = spectrum(words, 2 * k)
            if counts != spectrum(dual_words, 2 * k):
                continue
            if all(x.bit_count() % 4 == 0 for x in words):
                count_doubly_even_fsd += 1
                assert words == dual_words
                assert convex_order(counts, k)
            if min(x.bit_count() for x in words if x) >= 3:
                count_fsd_d3 += 1
                failures += not convex_order(counts, k)
        assert failures == 0
        print(f'STRENGTHENING k={k}: presentations={1 << (k*k)}; '
              f'FSD with d>=3={count_fsd_d3}; hinge failures={failures}; '
              f'doubly-even FSD={count_doubly_even_fsd} (all self-dual)')
    print('STRENGTHENING LIMIT: finite d>=3 checks do not settle all lengths; '
          'doubly-even linear FSD implies self-dual by polarization and dimension')


def main():
    fidelity_h8()
    family()
    witness()
    strengthenings()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
