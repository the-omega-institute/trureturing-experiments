#!/usr/bin/env python3
"""Exact finite checks of Gangloff--Nunez, arXiv:2605.04851v3, Question 4.20.

Subtraction is calculated by its defining meet, independently of the formula
for principal ideals used in the mathematical outline. No Lean proof is run.
"""
import itertools
import json


def subset(a, b):
    return a & b == a


def meet(values, top):
    out = top
    for value in values:
        out &= value
    return out


def ideals(n, strict):
    return [s for s in range(1 << n)
            if all(not (s >> j & 1) or (s >> i & 1) for i, j in strict)]


class IdealLattice:
    def __init__(self, n, strict):
        self.n = n
        self.strict = strict
        self.top = (1 << n) - 1
        self.elements = ideals(n, strict)
        self.maximal = {
            x: [m for m in self.elements if m != x and subset(m, x)
                and not any(m != u != x and subset(m, u) and subset(u, x)
                            for u in self.elements)]
            for x in self.elements}
        self.mu = {x: meet(ms, self.top) if ms else x
                   for x, ms in self.maximal.items()}
        self.subtraction = {}
        for x in self.elements:
            for z in self.elements:
                if subset(z, x):
                    eligible = [y for y in self.elements
                                if subset(y, x) and z | y == x]
                    assert eligible  # x itself is always eligible
                    result = meet(eligible, self.top)
                    assert result in self.elements
                    self.subtraction[x, z] = result
                    # Fidelity: source Lemma 4.2, z join (x-z) = x.
                    assert z | result == x
        for x, y in itertools.product(self.elements, repeat=2):
            assert x & y in self.elements and x | y in self.elements

    def boundary(self, x):
        sequence, strata, rank = [x], [], {}
        while self.mu[sequence[-1]] != sequence[-1]:
            stage = sequence[-1]
            residues = sorted({self.subtraction[stage, m]
                               for m in self.maximal[stage]})
            strata.append(residues)
            for s in residues:
                assert s not in rank  # fidelity: unique stratum, Notation 4.18
                rank[s] = len(strata) - 1
            sequence.append(self.mu[stage])
        assert self.maximal[sequence[-1]] == []
        for s, t in itertools.product(rank, repeat=2):
            # Fidelity: Proposition 4.19, strict reverse-rank monotonicity.
            if s != t and subset(s, t):
                assert rank[s] > rank[t]
        covers = [(s, t, rank[s], rank[t])
                  for s, t in itertools.product(sorted(rank), repeat=2)
                  if s != t and subset(t, s)
                  and not any(u != s and u != t and subset(t, u) and subset(u, s)
                              for u in rank)]
        return sequence, strata, rank, covers

    def all_subfamilies_check(self):
        """Check coframe identity and source's empty-filtered convention."""
        compact = set(self.elements)
        filtered_count = 0
        for flags in range(1 << len(self.elements)):
            family = [x for i, x in enumerate(self.elements) if flags >> i & 1]
            infimum = meet(family, self.top)
            assert infimum in self.elements
            for z in self.elements:
                assert z | infimum == meet((z | s for s in family), self.top)
            filtered = all(any(subset(c, a) and subset(c, b) for c in family)
                           for a, b in itertools.product(family, repeat=2))
            if filtered:
                filtered_count += 1
                if family:
                    assert infimum in family
                for z in self.elements:
                    if subset(infimum, z) and not any(subset(f, z) for f in family):
                        compact.discard(z)
        # The empty family witnesses failure of dual compactness exactly at top.
        assert compact == set(self.elements) - {self.top}
        # Every proper element has its singleton representation; top has empty.
        for z in self.elements:
            representation = [] if z == self.top else [z]
            assert all(f in compact for f in representation)
            assert meet(representation, self.top) == z
        return filtered_count, sorted(compact)


def all_small_posets():
    counts, lattices, boundary_cases = {}, 0, 0
    for n in range(5):
        edges = [(i, j) for i in range(n) for j in range(n) if i != j]
        count = 0
        for flags in range(1 << len(edges)):
            strict = {e for i, e in enumerate(edges) if flags >> i & 1}
            if any((j, i) in strict for i, j in strict):
                continue
            if any((i, k) not in strict for i, j in strict for jj, k in strict
                   if j == jj):
                continue
            count += 1
            lattice = IdealLattice(n, strict)
            for x in lattice.elements:
                lattice.boundary(x)
                boundary_cases += 1
            lattices += 1
        counts[n] = count
    assert counts == {0: 1, 1: 1, 2: 3, 3: 19, 4: 219}
    return dict(labelled_posets=counts, lattices=lattices, boundary_cases=boundary_cases)


def family(k):
    # 0=b, 1=a; 2,...,k+1=c_1,...,c_k.
    n = k + 2
    strict = {(0, 1)} | {(0, j) for j in range(2, n)}
    strict |= {(i, j) for i in range(2, n) for j in range(i + 1, n)}
    lattice = IdealLattice(n, strict)
    sequence, strata, rank, covers = lattice.boundary(lattice.top)
    chain = lambda j: 1 | sum(1 << i for i in range(2, j + 2))
    expected_ideals = sorted({0} | {chain(j) for j in range(k + 1)}
                             | {chain(j) | 2 for j in range(k + 1)})
    assert lattice.elements == expected_ideals
    assert len(lattice.elements) == 2 * k + 3
    assert sequence == [lattice.top] + [chain(k-i) for i in range(1, k+1)] + [0]
    assert strata == [[3, chain(k)]] + [[chain(k-i)] for i in range(1, k+1)]
    assert rank[3] == 0 and rank[1] == k
    assert (3, 1, 0, k) in covers
    assert k != 1 and any(rt != rs + 1 for _, _, rs, rt in covers)
    # Independent check of the principal-ideal residue formula, all ideals.
    for x in lattice.elements:
        maximal_points = [p for p in range(n) if x >> p & 1
                          and not any(i == p and x >> j & 1 for i, j in strict)]
        assert sorted(lattice.maximal[x]) == sorted(x ^ (1 << p) for p in maximal_points)
        for p in maximal_points:
            principal = (1 << p) | sum(1 << i for i, j in strict if j == p)
            assert lattice.subtraction[x, x ^ (1 << p)] == principal
    if k == 2:
        assert lattice.elements == [0, 1, 3, 5, 7, 13, 15]
        assert lattice.maximal == {0: [], 1: [0], 3: [1], 5: [1], 7: [3, 5],
                                   13: [5], 15: [7, 13]}
        assert [lattice.mu[x] for x in lattice.elements] == [0, 0, 1, 1, 1, 5, 5]
        assert covers == [(3, 1, 0, 2), (5, 1, 1, 2), (13, 5, 0, 1)]
        filtered, compact = lattice.all_subfamilies_check()
    else:
        filtered, compact = None, None
    return dict(k=k, ideal_count=len(lattice.elements), elements=lattice.elements,
                derivative=sequence + [sequence[-1]], residual_rank=len(strata),
                strata=strata, reverse_boundary_covers=covers,
                failing_cover=[3, 1], rank_gap=k,
                filtered_subfamilies_checked=filtered,
                dually_compact_allowing_empty_filtered=compact)


def main():
    print(json.dumps({'fidelity': 'Lemma 4.2 and Proposition 4.19',
                      'small_posets': all_small_posets()}, sort_keys=True))
    for k in range(2, 11):
        print(json.dumps(family(k), sort_keys=True))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
