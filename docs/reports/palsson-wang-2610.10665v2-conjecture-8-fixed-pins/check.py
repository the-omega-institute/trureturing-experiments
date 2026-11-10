#!/usr/bin/env python3
"""Exact finite checks for Palsson--Wang Conjecture 8 and a proposed shell family.

No finite experiment certifies an infinite dimension or an asymptotic density.
The accompanying preregistration supplies the unbounded mathematical argument.
Only integers, integer square roots and rational numbers are used as evidence.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import isqrt


def norm2(v):
    return sum(a * a for a in v)


def dist2(a, x):
    return sum((b - c) ** 2 for b, c in zip(a, x))


@lru_cache(None)
def pair_table(limit):
    ans = {}
    for c in range(isqrt(limit) + 1):
        for e in range(isqrt(limit - c * c) + 1):
            ans.setdefault(c * c + e * e, (c, e))
    return ans


@lru_cache(None)
def four_squares(n, limit):
    """Lexicographically first ordered nonnegative quadruple, exactly."""
    pairs = pair_table(limit)
    for a in range(isqrt(n) + 1):
        for b in range(isqrt(n - a * a) + 1):
            rem = n - a * a - b * b
            if rem in pairs:
                return (a, b) + pairs[rem]
    raise AssertionError((n, limit))


def fidelity():
    # Source Example 5.1, d=1: the pinned distance set of Z is the squares,
    # for every integer pin; this tests every truncated count for T <= 1000.
    cases = 0
    for x in range(-3, 4):
        for t in range(1001):
            h = isqrt(t)
            got = {dist2((a,), (x,)) for a in range(x - h, x + h + 1)}
            assert got == {k * k for k in range(h + 1)}
            assert len(got) == isqrt(t) + 1
            cases += 1
    # Source Example 5.1, d>=4: four squares produce every nonnegative value.
    for n in range(4097):
        w = four_squares(n, 4096)
        assert norm2(w) == n
    print(f'FIDELITY: Z pinned-square counts, {cases} cases; four-square coverage 0..4096')


def shell(d, j):
    m, s, radius = d - 4, 2 ** j, d * 2 ** j
    out = []
    for u in product(range(s, 2 * s), *([range(s)] * (m - 1))):
        assert norm2(u) < (d - 1) * s * s < radius * radius
        w = four_squares(radius * radius - norm2(u), radius * radius)
        a = u + w
        assert norm2(a) == radius * radius
        out.append(a)
    assert len(out) == s ** m
    assert len({a[:m] for a in out}) == s ** m
    return out


def finite_checks(d):
    m = d - 4
    js = [j for j in range(1, 4) if (2 ** j) ** m <= 40000]
    shells = {j: shell(d, j) for j in js}
    zero = (0,) * d
    basis = [zero] + [tuple(int(i == k) for i in range(d)) for k in range(d)]
    # The difference matrix of these affine points is literally the identity.
    assert [[basis[k + 1][i] - basis[0][i] for i in range(d)] for k in range(d)] == [list(x) for x in basis[1:]]
    assert len(set(basis)) == d + 1
    assert sum(Fraction(len(v), (2 ** j) ** m) for j, v in shells.items()) == len(js)
    # At alpha=m+1, the geometric scale factors sum to 1-2^-J.
    for J in range(1, 21):
        assert sum(Fraction(1, 2 ** j) for j in range(1, J + 1)) == 1 - Fraction(1, 2 ** J)
    # Small and remote cubes, including empty intersections. Square the mass
    # estimate to eliminate every Euclidean square root.
    cube_cases = 0
    for j, pts in shells.items():
        for a in [pts[0], pts[len(pts) // 2], pts[-1], zero]:
            for side in [1, 2, 3]:
                corner = tuple(v - side // 2 for v in a)
                r2 = d * side * side
                M2 = sum(max(c * c, (c + side) ** 2) for c in corner)
                hits = {k: sum(all(c <= v <= c + side for c, v in zip(corner, p)) for p in vs) for k, vs in shells.items()}
                if 4 * r2 < M2:
                    active = [k for k, n in hits.items() if n]
                    assert len(active) <= 2
                    for k in active:
                        assert hits[k] <= (side + 1) ** m
                        assert 4 * d * d * (2 ** k) ** 2 > M2
                    nu = sum(Fraction(n, (2 ** k) ** m) for k, n in hits.items())
                    C = 2 * (4 * d) ** m
                    alpha = m - 1 if m > 1 else 1
                    assert nu * nu * M2 ** alpha <= C * C * r2 ** alpha
                cube_cases += 1
    pins = basis + [tuple((-1) ** i * (i + 2) for i in range(d))]
    pins += [v[k] for v in shells.values() for k in [0, len(v) // 2, len(v) - 1]]
    distance_cases = 0
    for x in pins:
        X2 = norm2(x)
        B = isqrt(X2) + int(isqrt(X2) ** 2 != X2)
        for j, pts in shells.items():
            R = d * 2 ** j
            values = set()
            for a in pts:
                v = dist2(a, x)
                dot = sum(b * c for b, c in zip(a, x))
                assert v == R * R + X2 - 2 * dot
                assert dot * dot <= norm2(a) * X2
                assert R * R + X2 - 2 * B * R <= v <= R * R + X2 + 2 * B * R
                values.add(v)
            assert len(values) <= 4 * B * R + 1
            for T in [0, 1, R * R // 4, R * R, (R + B) ** 2]:
                if any(v <= T for v in values):
                    assert R <= isqrt(T) + B
            if x == zero:
                assert values == {R * R}
            distance_cases += len(pts)
    e1 = basis[1]
    for j, pts in shells.items():
        s, R = 2 ** j, d * 2 ** j
        values = {dist2(a, e1) for a in pts}
        assert values == {R * R + 1 - 2 * k for k in range(s, 2 * s)}
        assert len(values) == s and max(values) <= (R + 1) ** 2
        N = 2 * (R + 1) ** 2
        # |Q_N|=N+1 in ambient dimension one; a positive exact lower bound
        # for the SQUARE of the theta=1/2 density on the indicated subsequence.
        assert Fraction(s * s, N + 1) >= Fraction(1, 9 * d * d)
    theta = Fraction(3, 4) - Fraction(2, d)
    strict = Fraction(m, d) > Fraction(1, 2)
    assert strict == (d >= 9)
    if strict:
        assert Fraction(1, 2) < theta < Fraction(m, d)
        assert m > d * theta
    print(f'FAMILY d={d}: m={m}, shells={js}, points={sum(map(len, shells.values()))}, cube_cases={cube_cases}, distance_identities={distance_cases}, strict_refutation_range={strict}')


def main():
    fidelity()
    for d in range(5, 13):
        finite_checks(d)
    theta = Fraction(19, 36)
    assert 9 * theta == Fraction(19, 4) < 5
    assert theta > Fraction(1, 2)
    print('THRESHOLD d=9: theta=19/36, Hausdorff target=5 > 9*theta=19/4; theta>1/2')
    print('SCOPE: finite exact tests support the outline; they do not certify infinite dimensions or limits')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
