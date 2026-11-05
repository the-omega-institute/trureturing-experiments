#!/usr/bin/env python3
"""Exact finite CPS checks for Yang--Brandenburger, arXiv:2605.30017v2.

All arithmetic is integer or Fraction arithmetic. The search covers every
preorder on n <= 3 labelled states, every surjective ranking into 0,...,k-1
(k <= n), and the two positive weight vectors (1,...,1) and (1,...,n).
This is a restricted finite sample of real CPSs, not a proof of the claim.
"""
from fractions import Fraction as Q
from itertools import product
import json


def subsets(mask):
    sub = mask
    while True:
        yield sub
        if sub == 0:
            break
        sub = (sub - 1) & mask


def contained(a, b):
    return a & b == a


class CPS:
    def __init__(self, n, family, vectors):
        self.n = n
        self.full = (1 << n) - 1
        self.family = tuple(sorted(family))
        self.vectors = {g: tuple(Q(x) for x in vectors[g]) for g in self.family}
        self.table = {g: tuple(sum((v[x] for x in range(n) if e >> x & 1), Q(0))
                               for e in range(1 << n)) for g, v in self.vectors.items()}
        self.minimal = tuple(self.atom(x) for x in range(n))

    def atom(self, x):
        answer = self.full
        for g in self.family:
            if g >> x & 1:
                answer &= g
        return answer

    def p(self, g, e):
        return self.table[g][e]

    def certain(self, e):
        return sum(1 << x for x, m in enumerate(self.minimal) if self.p(m, e) == 1)

    def assessment(self, e, q):
        return sum(1 << x for x, m in enumerate(self.minimal) if self.p(m, e) == q)

    def check_axioms(self):
        assert self.family and 0 not in self.family
        assert all(any(g >> x & 1 for g in self.family) for x in range(self.n))
        for g in self.family:
            assert all(v >= 0 for v in self.vectors[g])
            assert sum(self.vectors[g]) == 1 and self.p(g, g) == 1
            for f in self.family:
                assert g | f in self.family
                assert not g & f or g & f in self.family
                if contained(f, g):
                    for e in subsets(f):
                        assert self.p(g, e) == self.p(g, f) * self.p(f, e)
        assert all(m in self.family and m >> x & 1 for x, m in enumerate(self.minimal))

    def one_closed(self):
        return all(self.p(g, f) != 1 or f in self.family
                   for g in self.family for f in subsets(g))

    def reflection(self):
        return all(self.p(m, self.assessment(e, self.p(m, e))) == 1
                   for e in range(1 << self.n) for m in self.minimal)

    def key(self):
        return self.family, tuple(self.vectors[g] for g in self.family)


def terminal(a, b, event, qa, qb):
    left, right = a.assessment(event, qa), b.assessment(event, qb)
    for _ in range(2 * a.n + 1):
        new = left & a.certain(right), right & b.certain(left)
        if new == (left, right):
            return left, right
        left, right = new
    raise AssertionError('Decreasing recursion did not stabilize')


def common_atom(a, b, x):
    answer = a.full
    for g in set(a.family) & set(b.family):
        if g >> x & 1:
            answer &= g
    return answer


def fidelity():
    # Source Example 1, Section 2: common certainty of disagreement without 1-closedness.
    a = CPS(4, (9, 6, 15), {9: (0, 0, 0, 1), 6: (0, 1, 0, 0), 15: (0, 0, 0, 1)})
    b = CPS(4, (7, 8, 15), {7: (0, 0, 1, 0), 8: (0, 0, 0, 1), 15: (0, 0, 0, 1)})
    a.check_axioms(); b.check_axioms()
    assert a.reflection() and b.reflection()
    assert all(a.vectors[common_atom(a, b, x)] == b.vectors[common_atom(a, b, x)]
               for x in range(4))
    assert a.assessment(2, Q(1)) == 6 and b.assessment(2, Q(0)) == 15
    left, right = terminal(a, b, 2, Q(1), Q(0))
    assert (left, right, left & right) == (6, 7, 6)
    assert not a.one_closed() and not b.one_closed()
    print('FIDELITY ' + json.dumps({'source': 'Example 1, Section 2',
          'A_terminal': ['b', 'c'], 'B_terminal': ['a', 'b', 'c'],
          'common_certainty': ['b', 'c'], 'q_A': '1', 'q_B': '0'}))

    # Source Example 2 has an independent, explicitly isolated hypothesis discrepancy.
    b2 = CPS(4, (3, 12, 15), {3: (Q(1, 2), Q(1, 2), 0, 0),
             12: (0, 0, Q(1, 2), Q(1, 2)), 15: (0, 0, Q(1, 2), Q(1, 2))})
    b2.check_axioms()
    assert b2.p(15, 13) == 1 and 13 not in b2.family and not b2.one_closed()
    print('SOURCE_BOUNDARY ' + json.dumps({'source': 'Example 2, Section 3',
          'finding': 'Displayed Bob family is not 1-closed',
          'conull_missing_event': ['a', 'c', 'd'], 'p_B_Omega': '1',
          'relevance': 'Not a counterexample to the target, which retains 1-closedness'}))


def nonvacuity():
    def uniform(g):
        return tuple(Q(int(bool(g >> x & 1)), g.bit_count()) for x in range(4))
    a = CPS(4, (3, 15), {g: uniform(g) for g in (3, 15)})
    b = CPS(4, (5, 15), {g: uniform(g) for g in (5, 15)})
    for c in (a, b):
        c.check_axioms()
        assert c.one_closed()
    assert terminal(a, b, 9, Q(1, 2), Q(1, 2)) == (15, 15)
    assert all(a.vectors[common_atom(a, b, x)] == b.vectors[common_atom(a, b, x)]
               for x in range(4))
    assert a.p(a.minimal[2], 3) == Q(1, 2)
    own = a.assessment(3, Q(1, 2))
    assert own == 12 and a.p(a.minimal[2], own) == Q(1, 2)
    assert not a.reflection() and not b.reflection()
    print('NONVACUITY ' + json.dumps({'states': 4, 'common_certainty': 'Omega',
          'q_A': '1/2', 'q_B': '1/2', 'reflection_A': False, 'reflection_B': False}))


def models(n):
    full = (1 << n) - 1
    unique = {}
    preorder_count = 0
    for bits in range(1 << (n * n)):
        row = tuple((bits >> (n * i)) & full for i in range(n))
        if not all(row[i] >> i & 1 for i in range(n)):
            continue
        if not all(not row[i] >> j & 1 or contained(row[j], row[i])
                   for i in range(n) for j in range(n)):
            continue
        preorder_count += 1
        family = tuple(g for g in range(1, full + 1)
                       if all(not g >> i & 1 or contained(row[i], g) for i in range(n)))
        for k in range(1, n + 1):
            for ranks in product(range(k), repeat=n):
                if set(ranks) != set(range(k)):
                    continue
                for weights in ((1,) * n, tuple(range(1, n + 1))):
                    vectors = {}
                    for g in family:
                        rank = min(ranks[x] for x in range(n) if g >> x & 1)
                        selected = tuple(x for x in range(n) if g >> x & 1 and ranks[x] == rank)
                        denominator = sum(weights[x] for x in selected)
                        assert denominator > 0
                        vectors[g] = tuple(Q(weights[x], denominator) if x in selected else Q(0)
                                           for x in range(n))
                    c = CPS(n, family, vectors)
                    if c.one_closed():
                        unique[c.key()] = c
    result = tuple(unique.values())
    for c in result:
        c.check_axioms()
    return result, preorder_count


def small_cases(n):
    collection, preorder_count = models(n)
    tests = common_cases = local_pairs = support_tests = averaging_tests = 0
    for c in collection:
        for event in range(1 << n):
            # Source Lemma 1: truth plus certainty implies knowledge.
            for x, m in enumerate(c.minimal):
                if event >> x & 1 and c.p(m, event) == 1:
                    assert contained(m, event)
            # Source Lemma 6: averaging constant posteriors over overlapping atoms.
            for g in c.family:
                values = {c.p(c.minimal[x], event) for x in range(n) if g >> x & 1}
                if len(values) == 1:
                    assert c.p(g, event) == next(iter(values))
                    averaging_tests += 1
    for a, b in product(collection, repeat=2):
        for w in range(n):
            h = common_atom(a, b, w)
            if a.vectors[h] != b.vectors[h]:
                continue
            local_pairs += 1
            support = sum(1 << y for y in range(n) if a.vectors[h][y] > 0)
            assert h == support | (1 << w)
            support_tests += 1
            for c in (a, b):
                for y in range(n):
                    if support >> y & 1:
                        m = c.minimal[y]
                        assert contained(m, h)
                        assert a.p(h, 1 << y) == a.p(h, m) * c.p(m, 1 << y)
                        assert c.p(m, 1 << y) > 0
            for event in range(1 << n):
                tests += 1
                qa, qb = a.p(a.minimal[w], event), b.p(b.minimal[w], event)
                left, right = terminal(a, b, event, qa, qb)
                if not (left & right) >> w & 1:
                    continue
                common_cases += 1
                joint = h & left & right
                for x in range(n):
                    if joint >> x & 1:
                        assert contained(a.minimal[x], joint)
                        assert contained(b.minimal[x], joint)
                assert contained(h, left & right)
                assert a.p(h, event) == qa == qb == b.p(h, event)
    return {'n': n, 'preorders': preorder_count, 'distinct_1_closed_models': len(collection),
            'locally_consistent_pair_state_tests': local_pairs, 'support_tests': support_tests,
            'event_state_tests': tests, 'common_certainty_cases': common_cases,
            'averaging_tests': averaging_tests, 'failures': 0}


def main():
    fidelity()
    nonvacuity()
    results = [small_cases(n) for n in (1, 2, 3)]
    for row in results:
        print('SMALL_CASES ' + json.dumps(row, sort_keys=True))
    print('TOTAL ' + json.dumps({'event_state_tests': sum(r['event_state_tests'] for r in results),
          'failures': 0, 'scope': 'All labelled preorders n<=3; surjective ranks; two weight choices'}))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
