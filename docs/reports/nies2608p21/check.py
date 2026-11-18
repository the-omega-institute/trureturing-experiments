#!/usr/bin/env python3
"""Exact finite checks for the proposed Nies Problem 2.1 counterexample mechanism.

These checks do not certify existence or omega-categoricity of the infinite
BEKP coding, and do not constitute a Lean proof of Problem 2.1.
"""
from fractions import Fraction
from itertools import permutations, product


def fidelity_rational_reversal():
    """Check the source's stated representative x -> -x on rational affine maps."""
    slopes = sorted({Fraction(p, q) for p in range(1, 5) for q in range(1, 5)})
    offsets = sorted({Fraction(p, q) for p in range(-4, 5) for q in range(1, 5)})
    maps = list(product(slopes, offsets))

    def compose(f, g):
        a, b = f
        c, d = g
        return a * c, a * d + b

    def reflected(f):
        a, b = f
        return a, -b

    points = [Fraction(p, q) for p in range(-4, 5) for q in range(1, 5)]
    for f in maps:
        a, b = f
        assert a > 0
        assert reflected(reflected(f)) == f
        for x in points:
            assert a * x - b == -(a * (-x) + b)
        for g in maps:
            assert reflected(compose(f, g)) == compose(reflected(f), reflected(g))
    translation = (Fraction(1), Fraction(1))
    assert reflected(translation) == (Fraction(1), Fraction(-1))
    assert reflected(translation) != translation
    # The two candidate normalizer actions on the pair-order types =,<,>.
    # Equality is fixed and the remaining two types can only be fixed or swapped.
    actions = [p for p in permutations(range(3)) if p[0] == 0]
    assert actions == [(0, 1, 2), (0, 2, 1)]
    print(f"FIDELITY rational reversal: {len(maps)} affine maps, "
          f"{len(maps) ** 2} composition pairs; two pair-order actions")


def synchronized_cover():
    """Exhaust all cover permutations, with E,F and an optional lifted order."""
    for n in range(1, 5):
        vertices = list(product(range(n), range(2)))
        pairs = list(product(range(2 * n), repeat=2))
        e = lambda x, y: vertices[x][0] == vertices[y][0]
        f = lambda x, y: vertices[x][1] == vertices[y][1]
        lt = lambda x, y: vertices[x][0] < vertices[y][0]
        actual, ordered = set(), set()
        for p in permutations(range(2 * n)):
            if not all(e(x, y) == e(p[x], p[y]) and f(x, y) == f(p[x], p[y])
                       for x, y in pairs):
                continue
            actual.add(p)
            if all(lt(x, y) == lt(p[x], p[y]) for x, y in pairs):
                ordered.add(p)
        expected = {tuple(2 * h[b] + (i ^ c) for b, i in vertices)
                    for h in permutations(range(n)) for c in range(2)}
        expected_ordered = {tuple(2 * b + (i ^ c) for b, i in vertices)
                            for c in range(2)}
        assert actual == expected
        assert ordered == expected_ordered
        print(f"COVER base_size={n}: Aut(E,F)={len(actual)}, "
              f"Aut(E,F,lifted_order)={len(ordered)}; exact product actions")


def abelian_shears():
    hom_pairs = 0
    for d in range(7):
        size = 1 << (d + 1)
        mask = (1 << d) - 1

        def shear(s, x):
            return x ^ (((s & (x & mask)).bit_count() & 1) << d)

        family = [tuple(shear(s, x) for x in range(size)) for s in range(1 << d)]
        assert len(set(family)) == 1 << d
        for s, alpha in enumerate(family):
            assert sorted(alpha) == list(range(size))
            assert all(alpha[alpha[x]] == x for x in range(size))
            for x, y in product(range(size), repeat=2):
                assert alpha[x ^ y] == (alpha[x] ^ alpha[y])
                hom_pairs += 1
            if s:
                j = (s & -s).bit_length() - 1
                assert shear(s, 1 << j) != 1 << j
                assert ((s & (1 << j)).bit_count() & 1) == 1
            for t in range(1 << d):
                assert all(alpha[family[t][x]] == family[s ^ t][x] for x in range(size))
        # This finite G is abelian, so its only inner automorphism is identity.
        assert all((g ^ x ^ g) == x for g, x in product(range(size), repeat=2))
        assert sum(a != tuple(range(size)) for a in family) == (1 << d) - 1
        print(f"SHEARS abelian d={d}: {len(family)} distinct outer classes; "
              "independence, homomorphism, involution and symmetric-difference law")
    print(f"SHEARS abelian homomorphism pairs checked={hom_pairs}")


def nonabelian_shears():
    s3 = list(permutations(range(3)))
    index = {p: i for i, p in enumerate(s3)}
    pmul = [[index[tuple(p[q[i]] for i in range(3))] for q in s3] for p in s3]
    identity = index[(0, 1, 2)]
    pinv = [next(j for j in range(6) if pmul[i][j] == identity and
                 pmul[j][i] == identity) for i in range(6)]
    assert any(pmul[a][b] != pmul[b][a] for a, b in product(range(6), repeat=2))
    for d in range(4):
        group = list(product(range(2), range(1 << d), range(6)))
        gi = {x: i for i, x in enumerate(group)}

        def mul(x, y):
            return x[0] ^ y[0], x[1] ^ y[1], pmul[x[2]][y[2]]

        def inv(x):
            return x[0], x[1], pinv[x[2]]

        def shear(s, x):
            return x[0] ^ ((s & x[1]).bit_count() & 1), x[1], x[2]

        inners = {tuple(gi[mul(mul(g, x), inv(g))] for x in group) for g in group}
        family = [tuple(gi[shear(s, x)] for x in group) for s in range(1 << d)]
        assert len(inners) == 6
        for s, a in enumerate(family):
            assert (a in inners) == (s == 0)
            for x, y in product(group, repeat=2):
                assert shear(s, mul(x, y)) == mul(shear(s, x), shear(s, y))
            for t, b in enumerate(family):
                composed = tuple(a[b[i]] for i in range(len(group)))
                assert composed == family[s ^ t]
                assert (composed in inners) == (s == t)
        print(f"SHEARS nonabelian d={d}: |G|={len(group)}, "
              f"|Inn(G)|={len(inners)}, outer shear classes={len(family)}")


def main():
    fidelity_rational_reversal()
    synchronized_cover()
    abelian_shears()
    nonabelian_shears()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
