#!/usr/bin/env python3
"""Exact checks for arXiv:2608.17990v2, Proposition 7.1 and Conjecture 7.2.

Python 3 and SymPy. All comparisons use integers, Fraction, or exact radicals.
This is evidence for a preregistration outline, not a Lean proof.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations
from math import lcm
import sympy as sp


def chirp(m):
    assert m >= 2 and m % 4 == 2
    a, b = Counter({m: 1}), Counter()
    for r in range(3, m - 2, 4):
        a[r] += 2
    for r in range(1, m, 4):
        b[r] += 2
    assert all(r > 0 and count > 0 for p in (a, b) for r, count in p.items())
    assert mass(a) == mass(b) == m * (m + 2) // 4
    return a, b


def mass(p):
    return sum(r * count for r, count in p.items())


def union(ms, padding=0):
    a, b = Counter(), Counter()
    for m in ms:
        c, d = chirp(m)
        a.update(c)
        b.update(d)
    a[1] += padding
    b[1] += padding
    assert padding >= 0 and mass(a) == mass(b)
    return a, b


def stairs(p, top):
    u, running = [0] * (top + 1), 0
    for s in range(top, 0, -1):
        running += s * p[s]
        u[s] = running
    return u


def profile(a, b):
    top = max(max(a), max(b))
    x, y = stairs(a, top), stairs(b, top)
    return [0] + [x[s] - y[s] for s in range(1, top + 1)]


def g(s):
    if s == 1:
        return 0
    a = s // 2
    return (-1) ** (a + 1) * 2 * a


def levels(n, w):
    assert n >= len(w) - 1 >= 2
    # Exact J = ceil(log_2 n), N = 2^J. Zero entries above support
    # contribute nothing; we still include every level through J.
    j = (n - 1).bit_length()
    e, blocks = [], []
    for ell in range(j + 1):
        sums = defaultdict(int)
        for s in range(1, len(w)):
            sums[(s - 1) // (1 << ell)] += w[s]
        assert all(0 <= b < (1 << j) // (1 << ell) for b in sums)
        e.append(sum(z * z for z in sums.values()))
        blocks.append(dict(sums))
    return e, blocks


def energy(e):
    # 2^(-ell/2) is 1/2^h for ell=2h, sqrt(2)/2^(h+1)
    # for ell=2h+1. Return the exact pair a,b for a+b*sqrt(2).
    a, b = Q(0), Q(0)
    for ell, z in enumerate(e):
        if ell % 2 == 0:
            a += Q(z, 1 << (ell // 2))
        else:
            b += Q(z, 1 << (ell // 2 + 1))
    return a, b


def sign2(a, b):
    """Exact sign of a+b*sqrt(2), for rational a,b."""
    if b == 0:
        return (a > 0) - (a < 0)
    if a == 0:
        return (b > 0) - (b < 0)
    if a > 0 and b > 0:
        return 1
    if a < 0 and b < 0:
        return -1
    d = a * a - 2 * b * b
    sign = (d > 0) - (d < 0)
    return sign if a > 0 else -sign


def ratio_below(n, norm, en, bound):
    a, b = en
    assert norm > 0 and n > 0 and a >= 0 and b >= 0 and bound > 0
    # Squaring is equivalent because both sides are positive.
    c = Q(n) * (a * a + 2 * b * b) - bound * bound * norm ** 4
    d = Q(n) * 2 * a * b
    return sign2(c, d) < 0


def check_nested(ms, padding=0):
    assert list(ms) == sorted(ms, reverse=True) and len(set(ms)) == len(ms)
    a, b = union(ms, padding)
    w = profile(a, b)
    assert w[1] == 0 and a != b
    assert w == [0] + [g(s) * sum(s <= m for m in ms) for s in range(1, len(w))]
    n = mass(a)
    assert n == sum(m * (m + 2) // 4 for m in ms) + padding
    norm = sum(abs(z) for z in w[1:])
    assert norm == sum(m * m // 2 for m in ms)
    e, blocks = levels(n, w)
    assert e[0] == sum((2 * j - 1) * (m ** 3 + 2 * m) // 3 for j, m in enumerate(ms, 1))
    assert e[1] == 2 * sum((2 * j - 1) * m for j, m in enumerate(ms, 1))
    assert all(z <= 4 * len(ms) ** 2 for z in e[2:])
    # Every individual layer has exactly one +2 block on each coarse level.
    for m in ms:
        c, d = chirp(m)
        _, bb = levels(n, profile(c, d))
        for row in bb[2:]:
            assert [z for z in row.values() if z != 0] == [2]
    assert all(sum(row.values()) == 2 * len(ms) for row in blocks[2:])
    en = energy(e)
    # Coarse tail <= (4+2*sqrt(2))*t^2, exact Q(sqrt(2)) comparison.
    ta, tb = energy([0, 0] + e[2:])
    t2 = len(ms) ** 2
    assert sign2(ta - 4 * t2, tb - 2 * t2) <= 0
    return n, norm, e, en


def main():
    print('Source: Liu-Wang-Tian-Wang, arXiv:2608.17990v2, Proposition 7.1 and Conjecture 7.2')
    print('Arithmetic: integers, Fraction, and SymPy exact algebra; SymPy', sp.__version__)
    fidelity = 0
    for m in range(2, 83, 4):
        n, norm, e, en = check_nested([m])
        a, b = chirp(m)
        w = profile(a, b)
        U = (m - 2) // 2
        assert e[0] == m * m + Q(4, 3) * U * (U + 1) * (2 * U + 1)
        assert e[1] == 2 * m and all(z == 4 for z in e[2:])
        assert sum(w) == 2
        tv = sum(abs(x - y) for x, y in zip(w[1:], w[2:] + [0]))
        assert tv == m * (m + 2) // 2 == 2 * n
        fidelity += 1
        if m == 10:
            assert (n, norm, e, en) == (30, 50, [340, 20, 4, 4, 4, 4], (Q(343), Q(23, 2)))
            print('Fidelity n=30: lambda=(10,7,7,3,3), mu=(9,9,5,5,1,1)')
            print('v=', w[1:], 'level_energies=', e, 'energy=343+(23/2)*sqrt(2)')
    print('Proposition 7.1 exact profiles, sizes, norms, variation and every level checked:', fidelity, 'chirps (m=2,6,...,82)')

    count = 0
    for t in (2, 3):
        for asc in combinations(range(2, 31, 4), t):
            check_nested(list(reversed(asc)))
            count += 1
    print('All nested pairs/triples from m in {2,6,...,30}:', count)

    for k in [2, 6, 10, 14, 18, 26, 50, 82, 126]:
        n, norm, e, en = check_nested([3 * k, k])
        assert n == (5 * k * k + 4 * k) // 2 and norm == 5 * k * k
        assert e[0] == 10 * k ** 3 + 4 * k and e[1] == 12 * k
        below = ratio_below(n, norm, en, Q(2, 3))
        if k >= 14:
            assert below
        if k >= 50:
            assert ratio_below(n, norm, en, Q(13, 20))
        print('Two layers: k=', k, 'n=', n, 'norm=', norm, 'E0=', e[0], 'E1=', e[1], 'energy=', en, 'ratio<2/3:', below)

    # A single rational bound establishes a uniform gap for all k>=50:
    # sqrt(n)<=8k/5, sqrt(2)<3/2, tail<28.
    assert Q(5, 2) + Q(2, 50) < Q(64, 25)
    assert Q(2) < Q(3, 2) ** 2
    gap_bound = Q(8, 125) * (10 + Q(13, 50 ** 2) + Q(28, 50 ** 3))
    assert gap_bound < Q(13, 20) < Q(2, 3)
    print('Uniform k>=50 upper bound:', gap_bound, '<13/20<2/3')

    k = sp.symbols('k', positive=True)
    two_majorant = (10 * k ** 3 + (4 + 6 * sp.sqrt(2)) * k + 16 + 8 * sp.sqrt(2)) * sp.sqrt((5 * k ** 2 + 4 * k) / 2) / (25 * k ** 4)
    two_minorant = (10 * k ** 3 + 4 * k) * sp.sqrt((5 * k ** 2 + 4 * k) / 2) / (25 * k ** 4)
    c = sp.sqrt(sp.Rational(2, 5))
    assert sp.simplify(sp.limit(two_majorant, k, sp.oo) - c) == 0
    assert sp.simplify(sp.limit(two_minorant, k, sp.oo) - c) == 0
    assert sp.Rational(2, 5) < sp.Rational(4, 9)
    m = sp.symbols('m', positive=True)
    n = m * (m + 2) / 4
    single_lo = ((m ** 3 + 2 * m) / 3) * sp.sqrt(n) / (m ** 4 / 4)
    single_hi = ((m ** 3 + 2 * m) / 3 + sp.sqrt(2) * m + 4 + 2 * sp.sqrt(2)) * sp.sqrt(n) / (m ** 4 / 4)
    assert sp.limit(single_lo, m, sp.oo) == sp.Rational(2, 3)
    assert sp.limit(single_hi, m, sp.oo) == sp.Rational(2, 3)
    print('Exact squeeze limits: source single chirp 2/3; two aligned chirps sqrt(2/5)<2/3')

    # Additional candidate screen: fixed-t reciprocal-odd layers, not needed
    # for the planned two-layer settlement. Padding is tested across a change of J.
    for t in range(1, 7):
        L = lcm(*range(1, 2 * t, 2))
        S = sum((Q(1, d * d) for d in range(1, 2 * t, 2)), Q(0))
        H = sum((Q(1, d) for d in range(1, 2 * t, 2)), Q(0))
        ms = [2 * L // d for d in range(1, 2 * t, 2)]
        n, norm, e, _ = check_nested(ms)
        assert n == L * L * S + L * H and norm == 2 * L * L * S
        assert e[0] == Q(8, 3) * L ** 3 * S + Q(4, 3) * t * L
        assert e[1] == 4 * t * L
        padding = (1 << (n - 1).bit_length()) + 1 - n
        check_nested(ms, padding)
        print('Extra fixed-t family: t=', t, 'L=', L, 'S_t=', S, 'n=', n, 'singleton padding=', padding)
    L, S, H, t = sp.symbols('L S H t', positive=True)
    many_hi = (sp.Rational(8, 3) * L ** 3 * S + (sp.Rational(4, 3) + 2 * sp.sqrt(2)) * t * L + (4 + 2 * sp.sqrt(2)) * t ** 2) * sp.sqrt(L ** 2 * S + L * H) / (4 * L ** 4 * S ** 2)
    assert sp.simplify(sp.limit(many_hi, L, sp.oo) - 2 / (3 * sp.sqrt(S))) == 0
    print('Extra fixed-t exact majorant limit: 2/(3*sqrt(S_t)); no claim of sharpness or true-limit existence')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
