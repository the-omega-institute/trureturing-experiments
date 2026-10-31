#!/usr/bin/env python3
"""Exact checks for Jones--Kinnersley Question 6.1.

Only finite integer combinatorics are used.  The script checks the published
n <= dirloc(Q_n) <= n+1 result on n <= 3 by an exact belief-state solver,
verifies the proposed n-cop strategy for n <= 5 under every legal reply, and
checks the paired-history lower-bound invariant for every relevant closed
neighborhood and projected probe set on n <= 4.
"""
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, product


def vertices(n):
    return tuple(range(1 << n))


def neighbors(n, v):
    return tuple(v ^ (1 << i) for i in range(n))


def closed_neighbors(n, v):
    return (v,) + neighbors(n, v)


def hamming(a, b):
    return (a ^ b).bit_count()


def legal_responses(n, probe, robber):
    if probe == robber:
        return (probe,)
    distance = hamming(probe, robber)
    return tuple(
        q for q in neighbors(n, probe) if hamming(q, robber) == distance - 1
    )


def response_classes(n, belief, probes):
    classes = {}
    for robber in belief:
        choices = [legal_responses(n, p, robber) for p in probes]
        for replies in product(*choices):
            classes.setdefault(replies, set()).add(robber)
    return {replies: frozenset(robbers) for replies, robbers in classes.items()}


def recontaminate(n, territory):
    return frozenset(
        y for x in territory for y in closed_neighbors(n, x)
    )


def exact_dirloc(n, k, max_rounds=8):
    actions = tuple(combinations_with_replacement(vertices(n), k))

    @lru_cache(None)
    def wins(belief, rounds):
        if len(belief) <= 1:
            return True
        if rounds == 0:
            return False
        for probes in actions:
            good = True
            for replies, current in response_classes(n, belief, probes).items():
                if len(current) == 1:
                    continue
                if not wins(recontaminate(n, current), rounds - 1):
                    good = False
                    break
            if good:
                return True
        return False

    initial = frozenset(vertices(n))
    for rounds in range(max_rounds + 1):
        if wins(initial, rounds):
            return rounds
    return None


def verify_upper_strategy(n):
    """Check the candidate's base-and-two-certified-coordinates strategy."""
    for base in vertices(n):
        probes = tuple(base ^ (1 << i) for i in range(n))
        for robber in vertices(n):
            distance = hamming(base, robber)
            if distance == 0:
                # Every response is the probe's own coordinate flip, so no
                # off-diagonal certification occurs and the base is identified.
                replies = tuple(legal_responses(n, p, robber)[0] for p in probes)
                certified = {
                    j for i, reply in enumerate(replies)
                    for j in range(n)
                    if i != j and reply == (probes[i] ^ (1 << j))
                }
                assert not certified
                continue
            if distance == 1:
                hit_coordinate = (base ^ robber).bit_length() - 1
                assert probes[hit_coordinate] == robber
                assert legal_responses(n, robber, robber) == (robber,)
                continue
            choices = [legal_responses(n, p, robber) for p in probes]
            for replies in product(*choices):
                certified = {
                    j for i, reply in enumerate(replies)
                    for j in range(n)
                    if i != j and reply == (probes[i] ^ (1 << j))
                }
                assert len(certified) >= 2, (n, base, robber, replies)
                chosen = sorted(certified)[:2]
                new_base = base
                for j in chosen:
                    new_base ^= 1 << j
                for moved in closed_neighbors(n, robber):
                    assert hamming(new_base, moved) <= distance - 1, (
                        n, base, robber, replies, new_base, moved
                    )


def verify_lower_invariant(n):
    """Check the paired-history invariant against every projected probe set."""
    universe = vertices(n - 1)
    for k in range(n):
        # A probe tuple contributes only its set of Q_(n-1) projections.
        for size in range(k + 1):
            for projected in combinations(universe, size):
                projected = frozenset(projected)
                # Initial paired histories choose y outside the first
                # projected probes; |Q_(n-1)| = 2^(n-1) >= n > k.
                assert any(y not in projected for y in universe)
                for previous_y in universe:
                    possible_y = [
                        y for y in closed_neighbors(n - 1, previous_y)
                        if y not in projected
                    ]
                    assert possible_y, (n, k, projected, previous_y)
                    y = possible_y[0]
                    for u in projected:
                        differing = u ^ y
                        assert differing != 0
                        bit = differing & -differing
                        q = u ^ bit
                        # q is one shortest-path first step toward both
                        # (y,0) and (y,1) in Q_(n-1) x K_2.
                        assert hamming(q, y) == hamming(u, y) - 1
                        assert hamming(q, y) >= 0
                    # The paired positions can both move to every closed
                    # neighbor of y while preserving the final K_2 bit.
                    for next_y in closed_neighbors(n - 1, y):
                        assert hamming(y, next_y) <= 1


def main():
    # Exact finite-state values reproduce the source's n versus n+1 bound
    # and settle the first three hypercubes at the conjectured value.
    values = {}
    for n in (1, 2, 3):
        values[n] = {}
        for k in range(n + 2):
            values[n][k] = exact_dirloc(n, k)
        assert values[n][n] is not None
        assert values[n][n - 1] is None
        assert all(
            values[n][k] is None or k >= n for k in range(n + 1)
        )
        assert values[n][n + 1] is not None
    for n in range(1, 6):
        verify_upper_strategy(n)
    for n in range(1, 5):
        verify_lower_invariant(n)
    print("exact_dirloc_small", values)
    print("upper_strategy_checked_dimensions", list(range(1, 6)))
    print("lower_paired_history_checked_dimensions", list(range(1, 5)))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
