#!/usr/bin/env python3
"""Exact certificates for Lin--Yoshida Conjecture 23; no numerical LP solver.

Fidelity: Proposition 10's three-dimensional example.
Settlement: certify a five-leaf witness, 110 generated binary trees, and
all coordinate partitions of sizes 3,...,7 having at least three classes.
Every number used as evidence is an integer or fractions.Fraction.
"""
from fractions import Fraction as Q
from itertools import combinations, product
from random import Random


def dtr(u, v):
    delta = [a - b for a, b in zip(u, v)]
    assert delta and len(u) == len(v)
    # Independent evaluation of source equations (1) and (2).
    value = max(delta) - min(delta)
    assert value == max(abs(a - b) for a in delta for b in delta)
    return value


def objective(x, sample):
    return sum((dtr(x, s) for s in sample), Q(0))


def canonical(x):
    return tuple(a - x[0] for a in x)


def certify(w, sample, labels, epsilon=Q(1)):
    """Check a balanced exact dual certificate and its equality rigidity.

    For row i select coordinates a_i in class i and b_i in class i+1.
    The affine lower bound r[a_i]-r[b_i] is bounded by the row range.
    Check that the sum has zero coefficient at EVERY ambient coordinate.
    Its constant is the objective at w, certifying global optimality.
    Tight row extrema give inequalities y_p <= y_q. A strongly connected
    graph of these inequalities forces every residual coordinate equal,
    certifying full ambient projective uniqueness, including competitors
    that split a replicated class. This is a finite symbolic certificate,
    not a grid search or a claim that a floating point LP found a singleton.
    """
    n, k = len(w), len(sample)
    assert k >= 3 and len(labels) == n and set(labels) == set(range(k))
    assert epsilon > 0
    representatives = [labels.index(i) for i in range(k)]
    coeff = [0] * n
    constant = Q(0)
    graph = [set([i]) for i in range(n)]
    for i, row in enumerate(sample):
        a, b = representatives[i], representatives[(i + 1) % k]
        residual = [w[e] - row[e] for e in range(n)]
        assert residual[a] == max(residual)
        assert residual[b] == min(residual)
        assert residual[a] - residual[b] == 2 * epsilon
        coeff[a] += 1
        coeff[b] -= 1
        constant += row[b] - row[a]
        for e in range(n):
            if residual[e] == residual[a]:
                graph[e].add(a)  # y[e] <= y[a]
            if residual[e] == residual[b]:
                graph[b].add(e)  # y[b] <= y[e]
    assert coeff == [0] * n
    assert constant == objective(w, sample) == 2 * k * epsilon
    for start in range(n):
        reached, pending = {start}, [start]
        while pending:
            for v in graph[pending.pop()]:
                if v not in reached:
                    reached.add(v)
                    pending.append(v)
        assert len(reached) == n
    assert len(set(map(canonical, sample))) == k
    deletion_gaps = []
    for i, s in enumerate(sample):
        others = sample[:i] + sample[i + 1:]
        gap = objective(s, others) - objective(w, others)
        assert gap == 2 * (k - 1) * epsilon > 0
        deletion_gaps.append(gap)
    assert all(dtr(a, b) == 4 * epsilon for a, b in combinations(sample, 2))
    # Quotient invariance, using independently varied row translations.
    translated = [tuple(x + (3 * i - 7) for x in s) for i, s in enumerate(sample)]
    assert objective(tuple(x + 11 for x in w), translated) == constant
    return constant, min(deletion_gaps)


def balanced(w, labels):
    return [tuple(x + (-1 if label == i else 1) for x, label in zip(w, labels))
            for i in range(max(labels) + 1)]


def ultrametric(v, n):
    pairs = list(combinations(range(n), 2))
    assert len(v) == len(pairs) and min(v) > 0
    distances = dict(zip(pairs, v))
    def d(i, j):
        return 0 if i == j else distances[tuple(sorted((i, j)))]
    # All ordered triples, including repeated indices and zero diagonals.
    assert all(d(i, k) <= max(d(i, j), d(j, k)) for i, j, k in product(range(n), repeat=3))
    # A resolved ultrametric determines a rooted triplet at every distinct triple.
    triplets = []
    for i, j, k in combinations(range(n), 3):
        options = [(d(i, j), (i, j, k)), (d(i, k), (i, k, j)), (d(j, k), (j, k, i))]
        least = min(x[0] for x in options)
        smaller = [x[1] for x in options if x[0] == least]
        assert len(smaller) == 1
        triplets.append(smaller[0])
    return tuple(triplets)


def partitions(n):
    """Restricted growth strings enumerate each set partition exactly once."""
    def extend(prefix, largest):
        if len(prefix) == n:
            if largest >= 2:
                yield prefix
            return
        for label in range(largest + 2):
            yield from extend(prefix + [label], max(largest, label))
    yield from extend([0], 0)


def generated_tree(n, rng):
    # Merge clusters; each internal node owns exactly its cross-child leaf pairs.
    clusters = [(frozenset([i]), 0, None) for i in range(n)]
    owner, levels, internal_edges = {}, [], []
    while len(clusters) > 1:
        i, j = sorted(rng.sample(range(len(clusters)), 2), reverse=True)
        left, right = clusters.pop(i), clusters.pop(j)
        idx = len(levels)
        height = max(left[1], right[1]) + 2 * rng.randint(3, 7)
        levels.append(height)
        for child in [left, right]:
            if child[2] is not None:
                internal_edges.append((idx, child[2]))
        for a in left[0]:
            for b in right[0]:
                owner[tuple(sorted((a, b)))] = idx
        clusters.append((left[0] | right[0], height, idx))
    assert len(levels) == n - 1
    assert min(levels) > 2
    assert all(levels[parent] - levels[child] > 2 for parent, child in internal_edges)
    labels = [owner[e] for e in combinations(range(n), 2)]
    w = tuple(levels[i] for i in labels)
    sample = balanced(w, labels)
    return w, sample, labels


def main():
    # Fidelity to Proposition 10, n=3, verbatim coordinates in the paper.
    fidelity = [(-1, 1, 1), (1, -1, 1), (1, 1, -1)]
    value, gap = certify((0, 0, 0), fidelity, [0, 1, 2])
    assert value == 6 and gap == 4
    print('FIDELITY Proposition 10 n=3: objective=6, unique class=[0], deletion gap=4')

    # Five-leaf, four-distinct-observation witness in order 12,13,14,15,23,24,25,34,35,45.
    w = (4, 20, 20, 20, 20, 20, 20, 12, 12, 6)
    sample = [(2, 20, 20, 20, 20, 20, 20, 12, 12, 6),
              (6, 20, 20, 20, 20, 20, 20, 14, 14, 8),
              (4, 20, 20, 20, 20, 20, 20, 10, 10, 6),
              (4, 20, 20, 20, 20, 20, 20, 12, 12, 4)]
    labels = [0, 1, 1, 1, 1, 1, 1, 2, 2, 3]
    value, gap = certify(w, sample, labels)
    topology = ultrametric(w, 5)
    assert all(ultrametric(s, 5) == topology for s in sample)
    assert all(max(s) == 20 for s in [w] + sample)
    assert len(set(w)) == 4 and canonical(w) != canonical((1,) * 10)
    assert objective(w, sample) == 8
    assert all(objective(s, sample[:i] + sample[i+1:]) == 12 for i, s in enumerate(sample))
    assert all(objective(w, sample[:i] + sample[i+1:]) == 6 for i in range(4))
    scaled_w = tuple(Q(x, 20) for x in w)
    scaled_s = [tuple(Q(x, 20) for x in s) for s in sample]
    assert certify(scaled_w, scaled_s, labels, Q(1, 20)) == (Q(2, 5), Q(3, 10))
    print('WITNESS N=5: objective=8, pair distance=4, deletion objectives=12 vs 6, gap=6')
    print('WITNESS: all ultrametric inequalities; identical resolved triplets; common root height=10')
    print('WITNESS: ambient uniqueness certificate; nonstar; normalization by 20 verified exactly')

    rng = Random(337)
    generated = 0
    for n in range(4, 15):
        for _ in range(10):
            center, rows, labels = generated_tree(n, rng)
            certify(center, rows, labels)
            topology = ultrametric(center, n)
            assert all(ultrametric(s, n) == topology for s in rows)
            # Restore a common root height by projectively translating each row.
            common_root = [tuple(x + max(center) - max(s) for x in s) for s in rows]
            certify(center, common_root, labels)
            assert all(ultrametric(s, n) == topology for s in common_root)
            assert all(max(s) == max(center) and min(s) > 0 for s in common_root)
            generated += 1
    print(f'GENERATED seed=337: {generated} binary trees, N=4,...,14, ten per N; all certificates passed')

    total = 0
    counts = []
    for n in range(3, 8):
        count = 0
        for labels in partitions(n):
            w = tuple(3 * e - e * e for e in range(n))
            certify(w, balanced(w, labels), labels)
            count += 1
        counts.append((n, count))
        total += count
    assert counts == [(3, 1), (4, 7), (5, 36), (6, 171), (7, 813)]
    print(f'PARTITIONS: {total} exact certificates; (coordinates, partitions)={counts}')
    print('BOUNDARY: finite exact certificates validate these cases; the universal theorem remains a Lean obligation')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
