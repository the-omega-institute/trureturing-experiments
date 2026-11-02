#!/usr/bin/env python3
"""Exact finite checks for Barman–Singh, arXiv:2606.21928v2.

Run with Python 3; standard library only. These finite checks do not prove the
unbounded conjecture. Cosets gH retain coefficient multiplicities.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, prod


def partitions(n, bound=None):
    if n == 0:
        yield ()
        return
    if bound is None:
        bound = n
    for a in range(min(n, bound), 0, -1):
        for rest in partitions(n - a, a):
            yield (a,) + rest


def cycle_type(p):
    seen = set()
    lengths = []
    for start in range(len(p)):
        if start in seen:
            continue
        a, length = start, 0
        while a not in seen:
            seen.add(a)
            length += 1
            a = p[a]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def class_size(shape):
    return factorial(sum(shape)) // prod(a ** m * factorial(m)
                                        for a, m in Counter(shape).items())


def coset_key(p, shape):
    """pH is determined by the images of the consecutive domain blocks."""
    blocks = []
    offset = 0
    for a in shape:
        blocks.append(tuple(sorted(p[offset:offset + a])))
        offset += a
    assert offset == len(p)
    return tuple(blocks)


def falling(a, j):
    return 0 if a < j else prod(range(a - j + 1, a + 1))


def evaluator_crosscheck():
    """Check the coset evaluator against literal pair multiplication."""
    total = 0
    for n in range(1, 5):
        group = list(permutations(range(n)))
        for shape in partitions(n):
            identity_key = coset_key(tuple(range(n)), shape)
            subgroup = [h for h in group if coset_key(h, shape) == identity_key]
            assert len(subgroup) == prod(factorial(a) for a in shape)
            for ctype in partitions(n):
                X = [p for p in group if cycle_type(p) == ctype]
                cosets = Counter(coset_key(x, shape) for x in X)
                literal = Counter(tuple(x[h[i]] for i in range(n))
                                  for x in X for h in subgroup)
                assert all(literal[g] == cosets[coset_key(g, shape)] for g in group)
                total += 1
    print(f"Evaluator: {total} class/subgroup cases, n=1..4; literal pairs agree")


def fidelity():
    """Source Example 3.8: the S9 classes (5,4) and (5,2,1,1)."""
    wanted = {(5, 4), (5, 2, 1, 1)}
    X = [p for p in permutations(range(9)) if cycle_type(p) in wanted]
    sizes = Counter(cycle_type(p) for p in X)
    assert sizes == Counter({(5, 4): 18144, (5, 2, 1, 1): 18144})
    assert all(sizes[c] == class_size(c) for c in wanted)
    for shape, expected_index, expected_r in [((7, 1, 1), 72, 504),
                                             ((6, 2, 1), 252, 144)]:
        histogram = Counter(coset_key(p, shape) for p in X)
        index = factorial(9) // prod(factorial(a) for a in shape)
        assert index == expected_index == len(histogram)
        assert set(histogram.values()) == {expected_r}
        assert len(X) * prod(factorial(a) for a in shape) == expected_r * factorial(9)
        print(f"Fidelity Example 3.8: lambda={shape}, cosets={index}, "
              f"every coefficient={expected_r}")
    fixed = Counter(sum(p[i] == i for i in range(9)) for p in X)
    moments = [sum(Fraction(w * falling(a, j), len(X)) for a, w in fixed.items())
               for j in range(3)]
    assert moments == [1, 1, 1]
    print(f"Fidelity moments j=0,1,2: {','.join(map(str, moments))}; |X|={len(X)}")


def small_conjecture():
    """Every allowed union, every lambda, n=5..7; no moment filter is used."""
    overall = 0
    for n in range(5, 8):
        classes = defaultdict(list)
        for p in permutations(range(n)):
            classes[cycle_type(p)].append(p)
        types = sorted(classes)
        sizes = [len(classes[c]) for c in types]
        assert all(sizes[i] == class_size(c) for i, c in enumerate(types))
        for shape in partitions(n):
            length = len(shape)
            if length < 5:
                continue
            cosets = {}
            rows = []
            for i, c in enumerate(types):
                for p in classes[c]:
                    key = coset_key(p, shape)
                    if key not in cosets:
                        cosets[key] = len(rows)
                        rows.append([0] * len(types))
                    rows[cosets[key]][i] += 1
            hsize = prod(factorial(a) for a in shape)
            index = factorial(n) // hsize
            assert len(rows) == index
            assert all(sum(row) == hsize for row in rows)
            checked = integral = 0
            for s in range(1, length - 1):
                for chosen in combinations(range(len(types)), s):
                    checked += 1
                    size = sum(sizes[i] for i in chosen)
                    r, remainder = divmod(size, index)
                    if remainder:
                        continue  # Uniform integer coefficients would force |X|=r*[G:H].
                    integral += 1
                    assert any(sum(row[i] for i in chosen) != r for row in rows), (
                        n, shape, [types[i] for i in chosen], r)
            overall += checked
            print(f"Conjecture: n={n}, lambda={shape}, unions={checked}, "
                  f"integral-candidate-r={integral}, tilings=0")
    assert overall == 8888
    print(f"Conjecture: all {overall} allowed union/lambda cases for n=5..7 pass; "
          "all 5<=k<=length(lambda) and all r>=0 are covered")


def dominance_and_support_checks():
    dominance = signs = 0
    for n in range(5, 26):
        for shape in partitions(n):
            if len(shape) < 5:
                continue
            for j in range(1, len(shape)):
                hook = (n - j,) + (1,) * j
                for a in range(1, max(len(hook), len(shape)) + 1):
                    assert sum(hook[:a]) >= sum(shape[:a])
                dominance += 1
    for a in range(51):
        for j in range(49):
            q = falling(a, j + 1) - falling(a, j)
            assert q == falling(a, j) * (a - j - 1)
            assert (q < 0) == (a == j)
            if a == j:
                assert q == -factorial(j)
            if a >= j + 2:
                assert q > 0
            signs += 1
    print(f"Hook dominance: {dominance} cases, 5<=n<=25, 1<=j<length(lambda)")
    print(f"Integer-support polynomial: {signs} exact sign/identity checks, "
          "0<=a<=50, 0<=j<=48")




def incidence_checks():
    """Check the Boolean-lattice commutator and transfer-matrix regularity."""
    identities = 0
    for n in range(2, 11):
        for r in range((n - 1) // 2 + 1):
            levels = list(combinations(range(n), r))
            for S_tuple in levels:
                S = frozenset(S_tuple)
                du = Counter()
                ud = Counter()
                for x in set(range(n)) - S:
                    T = S | {x}
                    for y in T:
                        du[tuple(sorted(T - {y}))] += 1
                for y in S:
                    R = S - {y}
                    for x in set(range(n)) - R:
                        ud[tuple(sorted(R | {x}))] += 1
                diff = {T: du[T] - ud[T] for T in du.keys() | ud.keys()
                        if du[T] != ud[T]}
                assert diff == {S_tuple: n - 2 * r}
                assert n - 2 * r > 0
                identities += 1
    transfers = 0
    for n in range(2, 8):
        for shape in partitions(n):
            if len(shape) < 2:
                continue
            # Move one point from the last block into the first (largest) block.
            a, b = shape[0], shape[-1]
            target = (a + 1,) + shape[1:-1] + ((b - 1,) if b > 1 else ())
            source_cosets = {coset_key(p, shape)
                             for p in permutations(range(n))}
            target_cosets = {coset_key(p, target)
                             for p in permutations(range(n))}
            rows = Counter()
            for P in source_cosets:
                for x in P[-1]:
                    first = tuple(sorted(P[0] + (x,)))
                    last = tuple(y for y in P[-1] if y != x)
                    Q = (first,) + P[1:-1] + ((last,) if last else ())
                    assert Q in target_cosets
                    rows[Q] += 1
            assert len(rows) == len(target_cosets)
            assert set(rows.values()) == {a + 1}
            assert len(source_cosets) * b == len(target_cosets) * (a + 1)
            transfers += 1
    print(f"Boolean-lattice commutator: {identities} exact basis-vector identities, "
          "2<=N<=10, 0<=r<N/2")
    print(f"One-box transfer: {transfers} Young-partition incidence matrices, "
          "degrees 2..7; row/column sums agree exactly")


if __name__ == '__main__':
    evaluator_crosscheck()
    fidelity()
    small_conjecture()
    dominance_and_support_checks()
    incidence_checks()
    print('ALL_CHECKS_OK')
