#!/usr/bin/env python3
"""Exact finite checks for KLM arXiv:2609.25132v1, Question 3.38.

Only integer arithmetic and finite sets are used. Finite checks are not a
proof of the universal topological assertions in the proposed refutation.
"""
from itertools import combinations, product


def fidelity_weak_degree():
    """Finite instances of the remark immediately after Definition 3.16.

    On X = G = Z/NZ with its regular action and basepoint 0, every
    {0}-definable Sigma-coloring has at most 2**len(Sigma) colors.
    Enumerate every Sigma of at most three clopens and every representation
    of length 1,...,4, including all repetitions of clopens.
    """
    count = 0
    for N in range(1, 5):
        points = tuple(range(N))
        clopens = tuple(frozenset(x for x in points if mask & (1 << x))
                        for mask in range(1 << N))
        for size in range(1, min(3, len(clopens)) + 1):
            for sigma in combinations(clopens, size):
                for n in range(1, 5):
                    for U in product(sigma, repeat=n):
                        # Definition 3.1, computed by actual translated sets.
                        colors = {tuple(int(0 in {(g + x) % N for x in u})
                                        for u in U) for g in points}
                        assert len(colors) <= 2 ** size
                        count += 1
    print(f'FIDELITY: post-Definition-3.16 weak finite degree; '
          f'N=1..4, |Sigma|=1..3, n=1..4; representations={count}; PASS')
    print('FIDELITY: Sigma empty has no positive-length representation; PASS')


def periodic_external_colors():
    """Compute translated clopen membership independently of the residue formula."""
    for r in range(2, 7):
        for m in range(2, 7):
            # y_i(t) = 1 iff t == i modulo m, and zero otherwise.
            points = [tuple(int(t == i) for t in range(m)) for i in range(m)]
            assert all(0 <= s < r for y in points for s in y)

            def shift(y, g):
                return tuple(y[(t - g) % m] for t in range(m))

            def c(g):
                # y in T^g U iff T^(-g)y in U, with U={x:x(0)=1}.
                return tuple(int(shift(y, -g)[0] == 1) for y in points)

            lo, hi = -3 * m, 3 * m
            cardinalities = []
            for h in range(lo, hi + 1):
                colors = {c(h + a) for a in range(m)}
                assert colors == set(points)
                assert len(colors) == m
                cardinalities.append(len(colors))
                for a in range(m):
                    assert c(h + a) == tuple(int((h + a) % m == i)
                                             for i in range(m))
            # All h follow from period m: this also enumerates every h mod m.
            for g in range(-3 * m, 3 * m + 1):
                assert c(g + m) == c(g)
            print(f'EXTERNAL: r={r}, m={m}, n={m}, Sigma={{x:x(0)=1}}, '
                  f'A=0..{m-1}, h={lo}..{hi}; '
                  f'min_colors={min(cardinalities)}, '
                  f'max_colors={max(cardinalities)}; PASS')
    # Negative control: distinct parameters are essential.
    for m in range(2, 7):
        repeated = {(int(g % m == 0),) * m for g in range(m)}
        assert len(repeated) == 2
    print('NEGATIVE CONTROL: repeating one periodic parameter gives only 2 colors; PASS')


def dense_word_and_orbit_checks():
    """Bounded checks of the explicit word enumeration and fixed-point mechanism."""
    total = 0
    for r in range(2, 5):
        words = [w for length in range(1, 5)
                 for w in product(range(r), repeat=length)]
        sequence = []
        locations = {}
        for w in words:
            locations[w] = len(sequence)
            sequence.extend(w)
        assert all(tuple(sequence[locations[w]:locations[w]+len(w)]) == w
                   for w in words)

        def x0(t):
            return sequence[t] if 0 <= t < len(sequence) else 0

        def T(g, t):
            return x0(t - g)

        # Every clopen supported in [-2,2] can be represented by an arbitrary
        # Boolean truth table. Making all its inputs zero checks all such tables
        # at once, including unions of cylinders and their complements.
        for parameters in [(0,), (-4, 0, 5), (3, 3, -2)]:
            for A in [(-3,), (-2, 0, 4), (0, 1, 2, 3)]:
                constraints = [t - b + a for b in parameters for a in A
                               for t in range(-2, 3)]
                h = -1 - max(constraints)
                assert all(T(b - h - a, t) == 0
                           for b in parameters for a in A for t in range(-2, 3))
                # Directly verify the source membership convention for these
                # orbit parameters and each one-coordinate cylinder.
                for b in parameters:
                    for a in A:
                        for t in range(-2, 3):
                            assert T(b, t + h + a) == T(b - h - a, t)
                total += 1
        print(f'DENSE WORDS: r={r}, all nonempty words of lengths 1..4 '
              f'occur; words={len(words)}, concatenation_length={len(sequence)}; PASS')
    print(f'ORBIT: simultaneous zero-window checks={total}, '
          'arbitrary clopen truth tables on [-2,2]; PASS')


def main():
    fidelity_weak_degree()
    periodic_external_colors()
    dense_word_and_orbit_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
