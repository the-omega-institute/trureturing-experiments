#!/usr/bin/env python3
"""Exact finite checks for Boichuk's trust-gossip absorption question.

All states and computations use integers. These checks are evidence for
the preregistered general argument, not a proof of the all-n assertion.
"""

from collections import deque
from itertools import permutations
from random import Random


def coordinates(n):
    return [(z, y) for z in range(n) for y in range(n) if z != y]


def decode(n, mask):
    rows = [0] * n
    for index, (z, y) in enumerate(coordinates(n)):
        if (mask >> index) & 1:
            rows[z] |= 1 << y
    return rows


def trusts(rows, z, y):
    assert z != y, "self-opinion is undefined"
    return (rows[z] >> y) & 1


def absorbing(rows):
    n = len(rows)
    return all(not trusts(rows, z, a) or
               trusts(rows, z, y) == trusts(rows, a, y)
               for a, z, y in permutations(range(n), 3))


def reachable(rows, remaining):
    reach = {z: {z} | {y for y in remaining if y != z and
                       trusts(rows, z, y)} for z in remaining}
    for a in sorted(remaining):
        for z in sorted(remaining):
            if a in reach[z]:
                reach[z] |= reach[a]
    return reach


def construct(initial):
    """Finish a sink component's entire ancestor basin before recursing."""
    rows = initial[:]
    remaining = set(range(len(rows)))
    trace = []

    def fire(a, z, y):
        assert len({a, z, y}) == 3
        assert trusts(rows, z, a) == 1
        assert trusts(rows, z, y) != trusts(rows, a, y)
        if trusts(rows, a, y):
            rows[z] |= 1 << y
        else:
            rows[z] &= ~(1 << y)
        trace.append((a, z, y))

    while remaining:
        original = rows[:]
        reach = reachable(original, remaining)
        root = min(remaining, key=lambda z: (len(reach[z]), z))
        core = reach[root]
        assert all(reach[x] == core for x in core)
        basin = {z for z in remaining if reach[z] & core}
        assert all(not trusts(original, z, y)
                   for z in remaining - basin for y in basin)

        # Distances and parents are fixed before any peripheral row is edited.
        distance = {z: 0 for z in core}
        frontier = deque(sorted(core))
        while frontier:
            w = frontier.popleft()
            for z in sorted(basin):
                if z not in distance and trusts(original, z, w):
                    distance[z] = distance[w] + 1
                    frontier.append(z)
        assert set(distance) == basin
        parent = {z: min(w for w in basin if w != z and
                         distance[w] == distance[z] - 1 and
                         trusts(original, z, w)) for z in basin - core}

        core_bits = sum(1 << y for y in core)
        for z in sorted(core):
            target_row = core_bits & ~(1 << z)
            while rows[z] & core_bits != target_row:
                # Find a two-edge shortcut with bit operations, rather than
                # rescanning every ordered triple after each addition.
                available = None
                for a in sorted(core - {z}):
                    missing = rows[a] & target_row & ~rows[z]
                    if trusts(rows, z, a) and missing:
                        y = (missing & -missing).bit_length() - 1
                        available = (a, z, y)
                        break
                assert available is not None
                fire(*available)
        assert all(trusts(rows, z, y) for z in core for y in core if z != y)
        expected_root = sum(1 << y for y in core if y != root)
        assert rows[root] == expected_root

        for z in sorted(basin - core, key=lambda x: (distance[x], x)):
            w = parent[z]
            assert rows[z] == original[z]
            assert trusts(rows, z, w)
            if not trusts(rows, z, root):
                assert w != root
                fire(w, z, root)
            for y in range(len(rows)):
                if y not in (z, root) and (trusts(rows, z, y) !=
                                          trusts(rows, root, y)):
                    fire(root, z, y)

        for z in basin:
            assert rows[z] == sum(1 << y for y in core if y != z)
        for z in remaining - basin:
            assert rows[z] == original[z]
            assert all(not trusts(rows, z, y) for y in basin)
        remaining -= basin
    return rows, trace


def check_construction(initial):
    final, trace = construct(initial)
    n = len(initial)
    # Replay independently in an off-diagonal dictionary representation.
    replay = {(z, y): (initial[z] >> y) & 1 for z, y in coordinates(n)}
    changed = set()
    for a, z, y in trace:
        assert len({a, z, y}) == 3
        assert replay[z, a] == 1
        assert replay[z, y] != replay[a, y]
        assert (z, y) not in changed
        changed.add((z, y))
        replay[z, y] = replay[a, y]
    assert all(replay[z, y] == ((final[z] >> y) & 1)
               for z, y in coordinates(n))
    assert all(replay[z, a] == 0 or replay[z, y] == replay[a, y]
               for a, z, y in permutations(range(n), 3))
    assert absorbing(final)
    assert len(trace) == len(changed) <= n * (n - 1)
    return len(trace)


def fidelity_flat_preparation(n):
    """Shortest flat path to absorption or an available progress move.

    This reproduces Section 10, Question 1 at n=3,4 by reverse BFS on
    the entire state graph. R counts reachable ordered DISTINCT pairs.
    """
    pairs = coordinates(n)
    index = {pair: i for i, pair in enumerate(pairs)}
    triples = [(index[z, a], index[z, y], index[a, y])
               for a, z, y in permutations(range(n), 3)]
    state_count = 1 << len(pairs)
    r_counts = []
    for mask in range(state_count):
        reach = reachable(decode(n, mask), set(range(n)))
        r_counts.append(sum(len(s) - 1 for s in reach.values()))
    reverse_flat = [[] for _ in range(state_count)]
    distance = [-1] * state_count
    frontier = deque()
    absorbing_count = 0
    for mask in range(state_count):
        successors = {mask ^ (1 << zy) for za, zy, ay in triples
                      if ((mask >> za) & 1) and
                      (((mask >> zy) ^ (mask >> ay)) & 1)}
        progress = False
        for after in successors:
            assert r_counts[after] <= r_counts[mask]
            if r_counts[after] < r_counts[mask]:
                progress = True
            else:
                reverse_flat[after].append(mask)
        if not successors:
            absorbing_count += 1
        if not successors or progress:
            distance[mask] = 0
            frontier.append(mask)
    while frontier:
        after = frontier.popleft()
        for before in reverse_flat[after]:
            if distance[before] == -1:
                distance[before] = distance[after] + 1
                frontier.append(before)
    assert -1 not in distance
    expected = {3: 1, 4: 2}[n]
    assert max(distance) == expected
    print(f"FIDELITY n={n} states={state_count} absorbing={absorbing_count} "
          f"max_shortest_flat_preparation={max(distance)} source={expected}", flush=True)


def sampled_states(n, rng):
    pairs = coordinates(n)
    yield [0] * n
    yield [sum(1 << y for y in range(n) if y != z) for z in range(n)]
    yield [1 << ((z + 1) % n) for z in range(n)]
    yield [(1 << (z + 1)) if z + 1 < n else 0 for z in range(n)]
    yield [sum(1 << y for y in range(1, n))] + [0] * (n - 1)
    yield [0] + [1] * (n - 1)
    for numerator in (1, 2, 4, 6, 8, 9):
        for _ in range(25):
            rows = [0] * n
            for z, y in pairs:
                if rng.randrange(10) < numerator:
                    rows[z] |= 1 << y
            yield rows


def main():
    for n in (3, 4):
        fidelity_flat_preparation(n)
    exhaustive_total = 0
    for n in range(1, 5):
        count = 1 << (n * (n - 1))
        maximum = max(check_construction(decode(n, mask)) for mask in range(count))
        exhaustive_total += count
        print(f"CONSTRUCTION exhaustive n={n} states={count} "
              f"max_constructed_length={maximum} bound={n * (n - 1)}", flush=True)
    rng = Random(260704023)
    sample_total = 0
    for n in (5, 6, 10, 20, 50):
        lengths = [check_construction(s) for s in sampled_states(n, rng)]
        sample_total += len(lengths)
        print(f"CONSTRUCTION sampled n={n} states={len(lengths)} "
              f"max_constructed_length={max(lengths)} bound={n * (n - 1)}", flush=True)
    print(f"TOTAL exhaustive={exhaustive_total} additional={sample_total} "
          "seed=260704023 exact_integer_arithmetic=True")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
