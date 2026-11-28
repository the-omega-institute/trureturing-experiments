"""Integral harmonics of 3C-I over all labelled three-generation spurion supports."""

from collections import Counter
from itertools import combinations_with_replacement, permutations, product

from sympy import Matrix, eye
from sympy.matrices.normalforms import hermite_normal_form


GENERATIONS = range(3)
WEINBERG = tuple(combinations_with_replacement(GENERATIONS, 2))
WEIGHT = (2, -1, -1)
FULL = (0, 0, 0)


def collapse(edges):
    """Canonical labels for left fields identified by column differences."""
    labels = list(GENERATIONS)
    for i, j in edges:
        left, right = sorted((labels[i], labels[j]))
        labels = [left if label == right else label for label in labels]
    return tuple(labels)


def yukawa_partition(mask):
    columns = [[i for i in GENERATIONS if mask & (1 << (3*i + j))] for j in GENERATIONS]
    return collapse((i, k) for rows in columns for i in rows for k in rows)


def quark_partition(up, down):
    return collapse((i, k) for i in GENERATIONS for k in GENERATIONS
                    if up[i] == up[k] or down[i] == down[k])


def lepton_order(labels, mask):
    """Zero denotes infinite order in the cokernel (no positive harmonic)."""
    adjacency = {label: set() for label in labels}
    for bit, (i, j) in enumerate(WEINBERG):
        if mask & (1 << bit):
            a, b = labels[i], labels[j]
            adjacency[a].add(b)
            adjacency[b].add(a)
    weights = {v: sum(w for i, w in enumerate(WEIGHT) if labels[i] == v)
               for v in adjacency}
    color = {}
    order = 1
    for start in adjacency:
        if start in color:
            continue
        color[start] = 1
        pending = [start]
        component = []
        bipartite = True
        while pending:
            v = pending.pop()
            component.append(v)
            for w in adjacency[v]:
                if w not in color:
                    color[w] = -color[v]
                    pending.append(w)
                elif color[w] == color[v]:
                    bipartite = False
        if bipartite and sum(color[v]*weights[v] for v in component):
            return 0
        if not bipartite and sum(weights[v] for v in component) % 2:
            order = 2
    return order


def raw_order(up, down, charged, weinberg):
    """Independent column HNF test in the original 16 field coordinates."""
    fields = eye(16)
    q, u, d, l, e = [[fields[:, offset+i] for i in GENERATIONS]
                     for offset in (0, 3, 6, 9, 12)]
    h = fields[:, 15]
    cols = []
    for mask, left, right, higgs in ((up, q, u, h), (down, q, d, -h), (charged, l, e, -h)):
        for i, j in product(GENERATIONS, repeat=2):
            if mask & (1 << (3*i+j)):
                cols.append(left[i]+right[j]+higgs)
    cols.extend(l[i]+l[j]+2*h for bit, (i, j) in enumerate(WEINBERG) if weinberg & (1 << bit))
    target = sum((w*(3*q[i]+l[i]) for i, w in enumerate(WEIGHT)), 0*h)
    basis = hermite_normal_form(Matrix.hstack(*cols) if cols else Matrix.zeros(16, 0))
    for k in (1, 2):
        if hermite_normal_form(basis.row_join(k*target)) == basis:
            return k
    assert basis.rank() < basis.row_join(target).rank()
    return 0


def minimum_support(invertible):
    """Minimum total nonzero entries, optionally requiring three perfect matchings."""
    sizes = {}
    for mask in range(1 << 9):
        if invertible and not any(all(mask & (1 << (3*i+j)) for i, j in enumerate(perm))
                                  for perm in permutations(GENERATIONS)):
            continue
        part = yukawa_partition(mask)
        sizes[part] = min(sizes.get(part, 10), mask.bit_count())
    quark_min = min(sizes[u]+sizes[d] for u, d in product(sizes, repeat=2)
                    if quark_partition(u, d) == FULL)
    minima = {}
    for part, size in sizes.items():
        for mask in range(1 << len(WEINBERG)):
            order = lepton_order(part, mask)
            total = quark_min + size + mask.bit_count()
            minima[order] = min(minima.get(order, 34), total)
    return minima[1], minima[2]


def census():
    representatives = {}
    counts = Counter()
    for mask in range(1 << 9):
        part = yukawa_partition(mask)
        representatives.setdefault(part, mask)
        counts[part] += 1
    quarks = Counter()
    quark_examples = {}
    for up, down in product(counts, repeat=2):
        part = quark_partition(up, down)
        quarks[part] += counts[up]*counts[down]
        quark_examples.setdefault(part, (representatives[up], representatives[down]))
    leptons = Counter()
    for part, charged in representatives.items():
        for w in range(1 << len(WEINBERG)):
            order = lepton_order(part, w)
            leptons[order] += counts[part]
            for qpart, (up, down) in quark_examples.items():
                assert raw_order(up, down, charged, w) == (order if qpart == FULL else 0)
    total = Counter({k: n*quarks[FULL] for k, n in leptons.items()})
    total[0] += (sum(quarks.values())-quarks[FULL])*sum(leptons.values())
    assert sum(total.values()) == 2**33
    assert total == {0: 2101520996, 1: 6197467374, 2: 290946222}
    assert minimum_support(False) == (5, 6)
    assert minimum_support(True) == (13, 14)
    print("Minimum primitive/doubled support: 5/6; with three perfect matchings: 13/14.")
    print("Yukawa partitions:", sorted(counts.items()))
    print("Connected quark supports:", quarks[FULL])
    print("Lepton supports by harmonic:", sorted(leptons.items()))
    print("All supports by harmonic:", sorted(total.items()))


if __name__ == "__main__":
    census()
