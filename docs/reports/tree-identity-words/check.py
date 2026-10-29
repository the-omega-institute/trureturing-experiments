#!/usr/bin/env python3
"""Exact finite checks of Shapiro, arXiv:2610.14373v1, Conjecture 1(iii).

No third-party packages. All counts, permutations and comparisons are integral.
The bounded check does not establish the conjecture for arbitrary trees.
"""
from collections import Counter, defaultdict
from itertools import product
from math import factorial


def require(test, message):
    if not test:
        raise AssertionError(message)


def prufer_edges(n, code):
    degree = [1] * n
    for v in code:
        degree[v] += 1
    edges = []
    for v in code:
        leaf = next(i for i in range(n) if degree[i] == 1)
        edges.append(tuple(sorted((leaf, v))))
        degree[leaf] -= 1
        degree[v] -= 1
    if n > 1:
        edges.append(tuple(i for i in range(n) if degree[i] == 1))
    return tuple(sorted(edges))


def tree_code(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)

    def rooted(v, parent):
        return '(' + ''.join(sorted(rooted(u, v) for u in adj[v] if u != parent)) + ')'

    return min(rooted(v, -1) for v in range(n))


def all_tree_shapes(n):
    if n == 1:
        return [()]
    representatives = {}
    for code in product(range(n), repeat=n - 2):
        edges = prufer_edges(n, code)
        representatives.setdefault(tree_code(n, edges), edges)
    return [representatives[k] for k in sorted(representatives)]


def append_swap(perm, edge):
    """Right multiplication: (sigma*tau)(x) = sigma(tau(x))."""
    q = list(perm)
    u, v = edge
    q[u], q[v] = q[v], q[u]
    return tuple(q)


def evaluate(n, edges, word):
    p = tuple(range(n))
    for a in word:
        p = append_swap(p, edges[a])
    return p


def inverse(p):
    q = [0] * len(p)
    for i, v in enumerate(p):
        q[v] = i
    return tuple(q)


def all_identity_onto_count(n, edges):
    """Meet in the middle counts ALL length-2e words using every edge.

    Half words are grouped by exact permutation and used-edge bitmask.
    A concatenation is an onto identity word iff the two permutations are
    inverse and the union of their masks is the full edge set. No occurrence
    multiplicity or noncrossing assumption enters this count.
    """
    e = len(edges)
    grouped = defaultdict(Counter)
    for half in product(range(e), repeat=e):
        mask = 0
        for a in half:
            mask |= 1 << a
        grouped[evaluate(n, edges, half)][mask] += 1
    total = 0
    for p, left in grouped.items():
        right = grouped.get(inverse(p), {})
        for ml, nl in left.items():
            for mr, nr in right.items():
                if ml | mr == (1 << e) - 1:
                    total += nl * nr
    return total


def local_noncrossing(edges, positions):
    for a, (i, j) in enumerate(positions):
        for b in range(a):
            k, l = positions[b]
            crossing = i < k < j < l or k < i < l < j
            if crossing and set(edges[a]) & set(edges[b]):
                return False
    return True


def check_double_words(n, edges):
    """Enumerate each word with exactly two occurrences per edge once."""
    e = len(edges)
    remaining = [2] * e
    positions = [[] for _ in edges]
    word = []
    counts = Counter()
    identity = tuple(range(n))

    def visit(p):
        if len(word) == 2 * e:
            lhs = p == identity
            rhs = local_noncrossing(edges, positions)
            require(lhs == rhs, ('criterion', n, edges, tuple(word), p, rhs))
            counts['words'] += 1
            counts['identity'] += lhs
            counts['local'] += rhs
            return
        for a in range(e):
            if remaining[a]:
                remaining[a] -= 1
                positions[a].append(len(word))
                word.append(a)
                visit(append_swap(p, edges[a]))
                word.pop()
                positions[a].pop()
                remaining[a] += 1

    visit(identity)
    require(counts['words'] == factorial(2 * e) // (2 ** e), ('multiset count', n, edges))
    return counts


def brute_all_onto(n, edges):
    e = len(edges)
    total = 0
    for word in product(range(e), repeat=2 * e):
        if len(set(word)) == e and evaluate(n, edges, word) == tuple(range(n)):
            total += 1
    return total


def source_name(n, edges):
    degree = Counter(v for edge in edges for v in edge)
    degs = sorted(degree.values(), reverse=True)
    e = len(edges)
    if e == 0:
        return 'K1'
    if e == 1:
        return 'K2'
    if max(degs) == 2:
        return 'P' + str(n)
    if max(degs) == e:
        return 'K1,' + str(e)
    if e == 4:
        return 'Ch'
    if degs[:2] == [4, 2]:
        return 'S4'
    if degs[:2] == [3, 3]:
        return 'DS'
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    center = next(v for v in range(n) if len(adj[v]) == 3)
    lengths = []
    for first in adj[center]:
        length, parent, cur = 1, center, first
        while len(adj[cur]) == 2:
            nxt = next(v for v in adj[cur] if v != parent)
            length, parent, cur = length + 1, cur, nxt
        lengths.append(length)
    return 'S(' + ','.join(map(str, sorted(lengths, reverse=True))) + ')'


def main():
    # Remark 8, the complete two tree tables preceding Conjecture 1.
    source = {'K2': 1, 'P3': 4, 'K1,3': 30, 'P4': 42,
              'K1,4': 336, 'Ch': 608, 'P5': 816,
              'K1,5': 5040, 'S4': 11220, 'DS': 14480,
              'S(3,1,1)': 18880, 'S(2,2,1)': 19680, 'P6': 25520}
    shape_counts = [1, 1, 1, 2, 3, 6]
    checked = 0
    seen_source = {}
    print('Source: Boris Shapiro, arXiv:2610.14373v1, Remark 8 and Conjecture 1(iii)')
    print('Convention: right multiplication, (sigma*tau)(x)=sigma(tau(x))')
    for n in range(1, 7):
        shapes = all_tree_shapes(n)
        require(len(shapes) == shape_counts[n - 1], ('tree shapes', n, len(shapes)))
        for edges in sorted(shapes, key=lambda es: source_name(n, es)):
            name = source_name(n, edges)
            total = all_identity_onto_count(n, edges)
            double = check_double_words(n, edges)
            require(total == double['identity'], ('identity word not double', n, edges))
            if n <= 5:
                require(brute_all_onto(n, edges) == total, ('brute/MITM disagreement', n, edges))
            if name in source:
                require(total == source[name], ('source table', name, total, source[name]))
                seen_source[name] = total
            print(f'e={n - 1} tree={name} edges={edges} '
                  f'double_words={double["words"]} identity_onto={total} '
                  f'local_noncrossing={double["local"]}')
            checked += double['words']
    require(set(seen_source) == set(source), 'source table coverage')
    print('Fidelity: all 13 published tree-table values reproduced exactly')
    print('Settlement check: all onto identity words are double and all double words obey the iff')
    print(f'Coverage: 14 tree shapes, 0..5 edges, {checked} double-occurrence words; '
          'all onto words covered by exact meet-in-the-middle counts')
    print('Cross-check: literal enumeration of every onto word for every tree with 0..4 edges')
    print('Boundary: finite verification only; no proof of the unbounded assertion')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
