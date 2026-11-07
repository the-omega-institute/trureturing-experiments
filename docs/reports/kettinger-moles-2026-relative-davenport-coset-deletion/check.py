#!/usr/bin/env python3
"""Exact finite checks of arXiv:2605.17595v2, Definitions 2.1 and 2.4.

A sequence is a list with repetitions. Only its full sum must lie in S.
A forbidden subsequence is any nonempty choice of positions, including all
positions. Computations below use integers modulo n and exact finite sets.
"""
from itertools import product
from math import gcd


def nonempty_sums(seq, n):
    """Independent literal enumeration of every nonempty positional subset."""
    return {
        sum(seq[i] for i in range(len(seq)) if mask & (1 << i)) % n
        for mask in range(1, 1 << len(seq))
    }


def exhaustive_profile(n):
    """Enumerate all zero-sum-free multisets, pruning only forbidden prefixes.

    Sorting is lossless in an abelian group: positional subsums and the full
    sum are permutation invariant. Every rejected prefix already contains a
    forbidden subsequence, hence all of its extensions are forbidden too.
    Every prefix at length n-1 is extended and checked; no longer valid
    branch is silently truncated.
    """
    lengths = [set() for _ in range(n)]
    witness = {}
    counts = [0] * (n + 1)
    nodes = 0

    def visit(seq, sums, first):
        nonlocal nodes
        nodes += 1
        assert len(seq) < n, (n, seq, "length bound failed")
        total = sum(seq) % n
        lengths[total].add(len(seq))
        witness[total, len(seq)] = seq
        counts[len(seq)] += 1
        for x in range(first, n):
            new_sums = sums | {x} | {(s + x) % n for s in sums}
            if 0 not in new_sums:
                visit(seq + (x,), new_sums, x)

    visit((), set(), 0)
    return lengths, witness, counts, nodes


def small_d(lengths, subset):
    return max((k for a in subset for k in lengths[a]), default=0)


def large_d(lengths, subset):
    # Minimal POSITIVE length with no zero-sum-free S-sum sequence.
    return next(k for k in range(1, len(lengths) + 1)
                if not any(k in lengths[a] for a in subset))


def main():
    profiles = {}
    for n in range(1, 13):
        lengths, witnesses, counts, nodes = exhaustive_profile(n)
        profiles[n] = lengths, witnesses
        actual = [small_d(lengths, {a}) for a in range(n)]
        expected = [n - gcd(n, a) for a in range(n)]
        assert actual == expected, (n, actual, expected)
        assert small_d(lengths, set(range(n))) == n - 1
        # Source Theorem 2.13: every subgroup coset in a cyclic group.
        cosets = 0
        for m in range(1, n + 1):
            if n % m:
                continue
            for g in range(1, m + 1):
                coset = {(g + j * m) % n for j in range(n // m)}
                assert small_d(lengths, coset) == n - gcd(g, m)
                cosets += 1
        # Lemma 2.5 and Theorem 2.6, all nonempty subsets for n <= 8.
        if n <= 8:
            for mask in range(1, 1 << n):
                subset = {a for a in range(n) if mask & (1 << a)}
                d = small_d(lengths, subset)
                assert (d == 0) == (subset == {0})
                assert large_d(lengths, subset) == d + 1
        # Surviving weaker form: H1 < H2, alpha=H1 and beta=H2.
        subgroup_pairs = 0
        for m2 in range(1, n + 1):
            if n % m2:
                continue
            h2 = set(range(0, n, m2))
            for m1 in range(m2 + 1, n + 1):
                if n % m1 or m1 % m2:
                    continue
                h1 = set(range(0, n, m1))
                assert h1 < h2
                assert small_d(lengths, h2) == small_d(lengths, h2 - h1)
                subgroup_pairs += 1
        print(f"FIDELITY C_{n}: singleton_profile={actual}; "
              f"cosets={cosets}; unshifted_subgroup_pairs={subgroup_pairs}; "
              f"zero_sum_free_multisets={nodes}; length_counts={counts}")

    # A second implementation exhausts ALL ordered C4 sequences of lengths
    # 0 through 4, using positional powersets directly rather than pruning.
    literal = [set() for _ in range(4)]
    tested = 0
    for length in range(5):
        for seq in product(range(4), repeat=length):
            tested += 1
            if 0 not in nonempty_sums(seq, 4):
                literal[sum(seq) % 4].add(length)
    assert literal == profiles[4][0]
    alpha = {1, 3}
    beta = set(range(4))
    retained = beta - alpha
    assert small_d(literal, beta) == 3
    assert small_d(literal, retained) == 2
    assert nonempty_sums((1, 1, 1), 4) == {1, 2, 3}
    assert nonempty_sums((1, 1), 4) == {1, 2}
    assert sum((1, 1)) % 4 in retained
    assert not set((1, 1)) <= retained  # Terms need not lie in S.
    print(f"C4_LITERAL: ordered_sequences={tested}; "
          f"lengths_by_total={[sorted(x) for x in literal]}; "
          "alpha=[1,3]; beta=[0,1,2,3]; retained=[0,2]; d_beta=3; d_retained=2")
    for n in range(4, 13, 2):
        lengths, witnesses = profiles[n]
        beta = set(range(n))
        alpha = set(range(1, n, 2))
        retained = beta - alpha
        assert len(alpha) >= 2 and len(retained) >= 2
        assert small_d(lengths, beta) == n - 1
        assert small_d(lengths, retained) == n - 2
        extremal_totals = {a for a in beta if n - 1 in lengths[a]}
        assert extremal_totals == {a for a in range(n) if gcd(a, n) == 1}
        assert extremal_totals <= alpha
        assert 0 not in nonempty_sums((1,) * (n - 2), n)
        assert sum((1,) * (n - 2)) % n in retained
        print(f"FAMILY C_{n}: d_beta={n-1}; d_even={n-2}; "
              f"extremal_totals={sorted(extremal_totals)}")
    # The Corollary 2.15 proof asserts this equality after adding p^b.
    p, a, b, g1 = 2, 2, 0, 1
    g1prime = g1 + p ** b
    left, right = gcd(g1prime, p ** b), gcd(g1prime, p ** a)
    assert (left, right) == (1, 2)
    print(f"COROLLARY_2_15_GCD_FAILURE: p={p}; a={a}; b={b}; "
          f"g1={g1}; g1prime={g1prime}; gcds={left},{right}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
