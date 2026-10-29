#!/usr/bin/env python3
"""Exact checks for Hakanen--Pavan arXiv:2610.35041v1, Conjecture 3.4.

No Lean proof is produced. None denotes unreachable distance; all arithmetic
is integral. The two solid-set predicates use independent quantifier forms.
"""
from collections import Counter, deque
from itertools import combinations, permutations, product

CERTIFICATES = {
    3: {2: [[[]]]},
    4: {4: [[[], [3]], [[2]]]},
    5: {
        8: [[[], [4], [3, 4]], [[3], [2, 3]], [[2], [2, 4]]],
        10: [[[], [5], [4, 5]], [[4], [3, 4]], [[3]]],
        12: [[[], [5], [3, 5]], [[4], [2, 4]], [[3], [3, 4]], [[2], [2, 5]]],
        40: [[[], [4], [3, 4]], [[3], [2, 3]], [[2], [2, 4]], [[1]]],
        41: [[[], [5], [1, 5]], [[3], [3, 5]], [[2], [2, 3]], [[1]]],
        76: [[[], [5], [3, 5]], [[4], [2, 4]], [[3], [1, 3]],
             [[2], [2, 5]], [[1], [1, 4]]],
    },
}
EXPECTED_ORBITS = {3: {2: 2}, 4: {4: 24},
                   5: {8: 120, 10: 120, 12: 120, 40: 40, 41: 120, 76: 24}}
EXPECTED_NONSTRONG = {3: 6, 4: 40, 5: 480}
EXPECTED_PROFILES = {(3, 2): 1, (4, 4): 3, (5, 8): 7, (5, 10): 6,
                     (5, 12): 9, (5, 40): 8, (5, 41): 8, (5, 76): 11}
LIMITS = {3: 1, 4: 2, 5: 5}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def adjacency(n, mask):
    a = [[False] * n for _ in range(n)]
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        a[i][j] = bool(mask & (1 << bit))
        a[j][i] = not a[i][j]
    return a


def distances(a):
    result = []
    for s in range(len(a)):
        ds = [None] * len(a)
        ds[s] = 0
        queue = deque([s])
        while queue:
            u = queue.popleft()
            for v in range(len(a)):
                if a[u][v] and ds[v] is None:
                    ds[v] = ds[u] + 1
                    queue.append(v)
        result.append(ds)
    return result


def lt(a, b):
    return a is not None and (b is None or a < b)


def solid_ordered(d, landmarks):
    return all(any(lt(d[s][u], d[s][v]) for s in landmarks)
               for u in range(len(d)) for v in range(len(d)) if u != v)


def solid_literal(d, landmarks):
    n = len(d)
    for v in range(n):
        for subset in range(1, 1 << n):
            if subset == 1 << v:
                continue
            witnessed = False
            for s in landmarks:
                finite = [d[s][x] for x in range(n)
                          if subset & (1 << x) and d[s][x] is not None]
                distance_to_set = min(finite) if finite else None
                if d[s][v] != distance_to_set:
                    witnessed = True
                    break
            if not witnessed:
                return False
    return True


def dimension(d):
    for k in range(len(d) + 1):
        for landmarks in combinations(range(len(d)), k):
            if solid_ordered(d, landmarks):
                return k
    raise AssertionError('The full vertex set must be solid')


def source_fidelity():
    # Theorem 2.1, checked from Definition 1.1 on every tournament and set n <= 4.
    comparisons = 0
    for n in range(1, 5):
        for mask in range(1 << (n * (n - 1) // 2)):
            d = distances(adjacency(n, mask))
            for subset in range(1 << n):
                landmarks = [s for s in range(n) if subset & (1 << s)]
                require(solid_literal(d, landmarks) == solid_ordered(d, landmarks),
                        ('Theorem 2.1 fidelity', n, mask, subset))
                comparisons += 1
    print('FIDELITY Theorem 2.1: all tournaments and all landmark sets, n=1..4;',
          comparisons, 'comparisons OK')
    # Proposition 3.2, exactly, on all labelled tournaments through five vertices.
    for n, lower in [(3, 3), (4, 3), (5, 4)]:
        histogram = Counter(dimension(distances(adjacency(n, mask)))
                            for mask in range(1 << (n * (n - 1) // 2)))
        require(min(histogram) >= lower, ('Proposition 3.2', n, histogram))
        print('FIDELITY Proposition 3.2:', n, 'dimension histogram', dict(sorted(histogram.items())), 'OK')
    # Proposition 3.3's transitive orientation has dimension n.
    for n in range(2, 12):
        d = distances(adjacency(n, 0))
        require(dimension(d) == n, ('Proposition 3.3', n))
    print('FIDELITY Proposition 3.3: transitive orientations n=2..11 OK')


def recode(a, permutation):
    return sum(int(a[permutation[i]][permutation[j]]) << bit
               for bit, (i, j) in enumerate(combinations(range(len(a)), 2)))


def terminal_component(d):
    # A vertex belongs to the last strong component iff it reaches no
    # vertex outside its own mutual-reachability class.
    return [i for i in range(len(d))
            if all(d[i][j] is None or d[j][i] is not None for j in range(len(d)))]


def check_chain_partition(chains, expected):
    flat = [tuple(x) for chain in chains for x in chain]
    require(len(flat) == len(set(flat)), 'duplicate profile')
    require(set(flat) == expected, 'profile coverage')
    for chain in chains:
        for a, b in zip(chain, chain[1:]):
            require(a != b and all(x <= y for x, y in zip(a, b)), 'not a strict chain')


def strong_certificates():
    for n, representatives in CERTIFICATES.items():
        orbit_sets = {m: {recode(adjacency(n, m), p) for p in permutations(range(n))}
                      for m in representatives}
        counts, nonstrong = Counter(), 0
        for mask in range(1 << (n * (n - 1) // 2)):
            d = distances(adjacency(n, mask))
            if any(None in row for row in d):
                nonstrong += 1
                continue
            hits = [m for m, orbit in orbit_sets.items() if mask in orbit]
            require(len(hits) == 1, ('classification', n, mask, hits))
            counts[hits[0]] += 1
        require(dict(counts) == EXPECTED_ORBITS[n], ('orbit sizes', n, counts))
        require(nonstrong == EXPECTED_NONSTRONG[n], ('nonstrong count', n))
        print('CLASSIFICATION', n, 'strong orbits', dict(sorted(counts.items())),
              'nonstrong', nonstrong, 'OK')
        for mask, encoded_chains in representatives.items():
            d = distances(adjacency(n, mask))
            radius = max(map(max, d))
            allowed = {x for x in product(range(1, 2 * radius), repeat=n)
                       if all(any(x[i] < d[i][j] for i in range(n)) for j in range(n))
                       and all(x[i] <= d[i][j] + x[j] for i in range(n) for j in range(n))}
            chains = [[tuple(2 if i + 1 in z else 1 for i in range(n)) for z in ch]
                      for ch in encoded_chains]
            check_chain_partition(chains, allowed)
            require(len(allowed) == EXPECTED_PROFILES[n, mask], ('profile count', n, mask))
            require(len(chains) <= LIMITS[n], ('chain bound', n, mask))
            print('STRONG PROFILE', n, mask, 'D=' + '/'.join(''.join(map(str, r)) for r in d),
                  'R', radius, 'profiles', len(allowed), 'chains', len(chains), 'OK')


def max_antichain(profiles):
    # These auxiliary nonstrong profile sets have at most eight members.
    best = 0
    for mask in range(1 << len(profiles)):
        selected = [x for i, x in enumerate(profiles) if mask & (1 << i)]
        if len(selected) <= best:
            continue
        if all(not all(a <= b for a, b in zip(x, y))
               and not all(b <= a for a, b in zip(x, y))
               for x, y in combinations(selected, 2)):
            best = len(selected)
    return best


def nonstrong_profiles():
    # For a terminal singleton B is impossible. Otherwise the proposed
    # terminal-component proof bounds every coordinate by 2; verify its
    # resulting relaxed profile space for every labelled nonstrong H.
    for n in [3, 4, 5]:
        cases, widths, terminal_sizes = 0, Counter(), Counter()
        for mask in range(1 << (n * (n - 1) // 2)):
            d = distances(adjacency(n, mask))
            if not any(None in r for r in d):
                continue
            terminal = terminal_component(d)
            require(len(terminal) in [1, 3, 4], ('terminal component size', n, mask))
            terminal_sizes[len(terminal)] += 1
            allowed = [x for x in product([1, 2], repeat=n)
                       if all(any(d[i][j] is None or x[i] < d[i][j]
                                  for i in range(n)) for j in range(n))
                       and all(d[i][j] is None or x[i] <= d[i][j] + x[j]
                               for i in range(n) for j in range(n))]
            width = max_antichain(allowed)
            require(width <= LIMITS[n], ('nonstrong width', n, mask, width))
            widths[width] += 1
            cases += 1
        print('NONSTRONG PROFILE', n, 'cases', cases, 'terminal sizes', dict(sorted(terminal_sizes.items())),
              'width histogram', dict(sorted(widths.items())), 'OK')


def direct_small_settlement():
    # A separate check on the original objects at the first conjectured order.
    histogram = Counter()
    for mask in range(1 << 15):
        d = distances(adjacency(6, mask))
        require(not any(solid_ordered(d, s) for s in combinations(range(6), 3)),
                ('Conjecture 3.4 n=6 counterexample', mask))
        histogram[dimension(d)] += 1
    print('DIRECT Conjecture 3.4 n=6: all 32768 labelled tournaments;',
          'dimension histogram', dict(sorted(histogram.items())), 'OK')


def main():
    source_fidelity()
    strong_certificates()
    nonstrong_profiles()
    direct_small_settlement()
    print('SCOPE: exact helper certificates for the general profile/terminal-component argument;')
    print('no enumeration of all tournaments at orders 7..11, no kernel proof claimed.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
