#!/usr/bin/env python3
"""Exact finite checks for Chen's Roman delta-complement equality question.

Python >= 3.10, standard library only. Read the six adjacent graph6 datasets;
see README.md for their provenance. No floating point and no Lean proofs.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

DATA = {
    3: (4, 'aefbaa12a956ed1f415fa897c455185134275a89a57ce1ef7d38f771c0d9129e'),
    4: (11, 'c98e9d5ed38843ff55effb40b3ccc33e348027626268109e4811b64f4d7455ff'),
    5: (34, '9481bb6227eb9eb8e3de2a9d98ba29428dce99fb640bdb8f5752914b40ddc0fb'),
    6: (156, '3e854c54447bbde07316b94338a0b39ddb7d2c094dd63fb8867df829fdcc8043'),
    7: (1044, '5a878f1a46606da4733df250bf9a412d6ae0b3c73c2112e25d66ba00ae565b56'),
    8: (12346, '546a249902101c97d3aa590f93e53366854bd0a6f405aa59bdb32d25c57f845a'),
}


def graph(n, edges):
    a = [0] * n
    for u, v in edges:
        assert 0 <= u < n and 0 <= v < n and u != v
        a[u] |= 1 << v
        a[v] |= 1 << u
    return tuple(a)


def graph6(line):
    b = line.strip()
    n = b[0] - 63
    assert 0 <= n <= 8
    bits = [(x - 63 >> j) & 1 for x in b[1:] for j in range(5, -1, -1)]
    pairs = [(i, j) for j in range(1, n) for i in range(j)]
    assert len(bits) >= len(pairs) and not any(bits[len(pairs):])
    return graph(n, (e for e, bit in zip(pairs, bits) if bit))


def complement(a):
    full = (1 << len(a)) - 1
    return tuple(full ^ row ^ (1 << i) for i, row in enumerate(a))


def degree_classes(a):
    classes = {}
    for i, row in enumerate(a):
        classes.setdefault(row.bit_count(), []).append(i)
    return list(classes.values())


def delta(a):
    deg = [row.bit_count() for row in a]
    return graph(len(a), ((u, v) for u, v in combinations(range(len(a)), 2)
                         if bool(a[u] & (1 << v)) ^ (deg[u] == deg[v])))


def induced(a, vertices):
    return tuple(sum(1 << j for j, v in enumerate(vertices) if a[u] & (1 << v))
                 for u in vertices)


def gamma(a):
    # For a fixed set S of 2-labelled vertices, the cheapest completion puts
    # 1 exactly outside N[S]. Enumerate every S, including the empty set.
    n = len(a)
    closed = [0] * (1 << n)
    best = n
    for s in range(1, 1 << n):
        bit = s & -s
        v = bit.bit_length() - 1
        closed[s] = closed[s ^ bit] | bit | a[v]
        best = min(best, 2 * s.bit_count() + n - closed[s].bit_count())
    return best


def is_roman(a, f):
    return all(f[v] != 0 or any(f[u] == 2 and (a[v] >> u) & 1
                              for u in range(len(a))) for v in range(len(a)))


def owner_k_one(a, f):
    return all(f[v] != 0 or
               1 + sum(1 for u in range(len(a)) if (a[v] >> u) & 1 and f[u] > 0)
               <= sum(f[u] for u in range(len(a)) if (a[v] >> u) & 1)
               for v in range(len(a)))


def gamma_ternary(a, check_bridge=False):
    best = len(a)
    for f in product(range(3), repeat=len(a)):
        valid = is_roman(a, f)
        if check_bridge:
            assert valid == owner_k_one(a, f), (a, f)
        if valid:
            best = min(best, sum(f))
    return best


def components(a):
    todo = set(range(len(a)))
    answer = []
    while todo:
        seen = {min(todo)}
        stack = list(seen)
        while stack:
            u = stack.pop()
            for v in list(todo - seen):
                if (a[u] >> v) & 1:
                    seen.add(v)
                    stack.append(v)
        todo -= seen
        answer.append(sorted(seen))
    return answer


def is_cycle(a, n):
    return len(a) == n and all(x.bit_count() == 2 for x in a) and len(components(a)) == 1


def is_matching_graph(a):
    return all(row.bit_count() == 1 for row in a)


def local_extremal(a):
    return (len(a) <= 2 or is_cycle(a, 5) or is_matching_graph(a)
            or is_matching_graph(complement(a)))


def structural(a):
    deg = [row.bit_count() for row in a]
    return (all(not ((a[u] >> v) & 1) or deg[u] == deg[v]
                for u, v in combinations(range(len(a)), 2))
            and all(local_extremal(induced(a, c)) for c in degree_classes(a)))


def explicit_family(a):
    z = 0
    cycle_choice = 0
    cp_sizes = set()
    for vertices in components(a):
        c = induced(a, vertices)
        if len(c) == 1:
            z += 1
        elif len(c) == 2 and is_matching_graph(c):
            pass
        elif is_cycle(c, 4) or is_cycle(c, 5):
            cycle_choice += 1
        elif len(c) >= 6 and is_matching_graph(complement(c)):
            if len(c) in cp_sizes:
                return False
            cp_sizes.add(len(c))
        else:
            return False
    return bool(a) and z <= 2 and cycle_choice <= 1


def bounds(a):
    sizes = [len(c) for c in degree_classes(a)]
    n, k = len(a), len(sizes)
    return n + 3 * k - 2 * sizes.count(1) - sizes.count(2)


def union(*graphs):
    a = []
    shift = 0
    for g in graphs:
        a.extend(row << shift for row in g)
        shift += len(g)
    return tuple(a)


def main():
    edge = graph(2, [(0, 1)])
    c4 = graph(4, [(i, (i + 1) % 4) for i in range(4)])
    example = union(edge, edge, c4)
    assert (gamma(example), gamma(delta(example)), bounds(example)) == (7, 7, 14)
    print('FIDELITY Example 2: gamma_R(2K2 disjoint C4)=7; gamma_R(delta)=7; refined bound=14')

    # Enumeration independent of the downloaded isomorphism-class datasets.
    labelled = 0
    for n in range(1, 6):
        pairs = list(combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            a = graph(n, (e for j, e in enumerate(pairs) if (mask >> j) & 1))
            assert gamma(a) == gamma_ternary(a)
            eq = gamma(a) + gamma(delta(a)) == bounds(a)
            assert eq == structural(a) == explicit_family(a), (n, mask)
            labelled += 1
    print(f'INDEPENDENT labelled graphs n=1..5: {labelled}; ternary and subset minima agree; classification agrees')

    root = Path(__file__).resolve().parent
    total = sharp_total = 0
    slack_sum = 0
    bridge_count = 0
    for n, (expected_count, digest) in DATA.items():
        data = (root / f'graph{n}.g6').read_bytes()
        assert sha256(data).hexdigest() == digest
        lines = data.splitlines()
        assert len(lines) == expected_count and len(set(lines)) == expected_count
        sharp = 0
        for line in lines:
            a = graph6(line)
            assert len(a) == n and complement(complement(a)) == a
            classes = degree_classes(a)
            ad = delta(a)
            for c in classes:
                assert induced(ad, c) == complement(induced(a, c))
            g, gd, b = gamma(a), gamma(ad), bounds(a)
            assert g + gd <= b, line
            eq = g + gd == b
            assert eq == structural(a) == explicit_family(a), line
            assert local_extremal(a) == (g + gamma(complement(a)) == n + 3), line
            if n <= 5:
                assert g == gamma_ternary(a, check_bridge=True)
                bridge_count += 3 ** n
            if n <= 6:
                assert gd == gamma_ternary(ad)
            slack_sum += b - g - gd
            sharp += eq
        assert sharp == {3: 1, 4: 3, 5: 3, 6: 6, 7: 5, 8: 8}[n]
        print(f'UNLABELLED n={n}: graphs={len(lines)}; bound_violations=0; equality={sharp}; classification_mismatches=0')
        total += len(lines)
        sharp_total += sharp
    assert (total, sharp_total) == (13595, 26)
    assert (100 * slack_sum + total // 2) // total == 791
    print(f'FIDELITY Table 7: total_graphs={total}; equality_count={sharp_total}; counts=(1,3,3,6,5,8)')
    print(f'FIDELITY Table 6: exact mean refined slack={Fraction(slack_sum, total)}; rounds to 7.91')
    print(f'BRIDGE finite labels checked={bridge_count}; owner IsRomanDom at k=1 agrees with ordinary Roman predicate')
    print('SCOPE: finite diagnostics only; the universal classification and bridge remain unproved in Lean')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
