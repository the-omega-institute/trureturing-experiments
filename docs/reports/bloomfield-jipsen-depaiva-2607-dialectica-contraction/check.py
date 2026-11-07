#!/usr/bin/env python3
"""Exact finite checks for the Dialectica contraction preregistration.

Only Python's standard library is used. Finite checks do not certify the
infinite family or the categorical universal properties.
"""

from itertools import product
from math import prod


def dle(a, b):
    ar, u, x = a
    br, v, y = b
    return u <= v and min(u, y) <= x and min(ar, u, y) <= min(br, u, y)


def tensor(a, b):
    return tuple(min(x, y) for x, y in zip(a, b))


def bound(objects, a, b, upper):
    candidates = [c for c in objects if
                  (dle(a, c) and dle(b, c) if upper else dle(c, a) and dle(c, b))]
    best = [c for c in candidates if all(
        dle(c, d) if upper else dle(d, c) for d in candidates)]
    assert len(best) == 1
    return best[0]


def fidelity():
    # Section 2, displayed tensor table for D'(2), in the source's order.
    objects = [(0, 0, 0), (0, 1, 1), (1, 1, 1), (0, 1, 0)]
    expected = [[0, 0, 0, 0], [0, 1, 1, 3], [0, 1, 2, 3], [0, 3, 3, 3]]
    actual = [[objects.index(tensor(a, b)) for b in objects] for a in objects]
    assert actual == expected
    assert all(tensor(a, a) == a for a in objects)
    # The source also states that D'(3) has ten elements and is nondistributive.
    three = [(a, u, x) for a, u, x in product(range(3), repeat=3) if a <= x <= u]
    assert len(three) == 10
    join = lambda a, b: bound(three, a, b, True)
    meet = lambda a, b: bound(three, a, b, False)
    failures = [(a, b, c) for a, b, c in product(three, repeat=3)
                if meet(a, join(b, c)) != join(meet(a, b), meet(a, c))]
    assert failures
    print('FIDELITY: Dprime(2) tensor table =', actual)
    print('FIDELITY: Dprime(3) cardinality = 10; distributivity failures =', len(failures))
    print('FIDELITY: first distributivity witness =', failures[0])
    # Section 4's membership-dependent Set contraction selector, including empties.
    relation_count = 0
    for usize in range(3):
        for xsize in range(4):
            pairs = list(product(range(usize), range(xsize)))
            for bits in product((False, True), repeat=len(pairs)):
                relation = {p for p, bit in zip(pairs, bits) if bit}
                for u, x, y in product(range(usize), range(xsize), range(xsize)):
                    f = x if (u, x) not in relation else y
                    assert (u, f) not in relation or (
                        (u, x) in relation and (u, y) in relation)
                relation_count += 1
    print('FIDELITY: source Set selector checked on', relation_count, 'relations')


def allowed(xset, aset, pair):
    x, y = pair
    return {z for z in xset if z not in aset or (x in aset and y in aset)}


def obstruction(p, q, r):
    pset = set(range(p))
    qset = set(range(p, p + q))
    rset = set(range(p + q, p + q + r))
    xset = pset | qset | rset
    stages = [rset, pset | rset, qset | rset]
    assert all(aset and aset < xset for aset in stages)
    assert stages[0] <= stages[1] and stages[0] <= stages[2]
    pairs = list(product(sorted(xset), repeat=2))
    # Identity structure maps make every natural transformation a shared function.
    choices = [set.intersection(*(allowed(xset, a, pair) for a in stages)) for pair in pairs]
    natural_selector_count = prod(map(len, choices))
    assert natural_selector_count == 0
    pair = (min(pset), min(qset))
    assert allowed(xset, stages[1], pair) == qset
    assert allowed(xset, stages[2], pair) == pset
    assert not pset & qset
    # Each stage separately admits selectors. Compatibility is the obstruction.
    stage_counts = [prod(len(allowed(xset, a, pair)) for pair in pairs) for a in stages]
    assert all(n > 0 for n in stage_counts)
    return natural_selector_count, stage_counts


def exact_witness():
    xset = set(range(3))
    stages = [{2}, {0, 2}, {1, 2}]
    pairs = list(product(range(3), repeat=2))
    count = 0
    stage_counts = [0, 0, 0]
    # Exhaust all 3^9 shared functions independently of the product-of-choices test.
    for values in product(range(3), repeat=9):
        good = [all(z in allowed(xset, a, pair) for pair, z in zip(pairs, values))
                for a in stages]
        stage_counts = [n + int(ok) for n, ok in zip(stage_counts, good)]
        count += int(all(good))
    assert count == 0
    assert obstruction(1, 1, 1) == (count, stage_counts)
    print('EXACT WITNESS: exhaustively checked 19683 natural candidate functions; accepted =', count)
    print('EXACT WITNESS: separately admissible stage-function counts =', stage_counts)
    print('EXACT WITNESS: at (0,1), leaf s allows {1}; leaf t allows {0}')


def family():
    triples = list(product(range(1, 5), repeat=3))
    for sizes in triples:
        assert obstruction(*sizes)[0] == 0
    print('FAMILY: all', len(triples), 'positive block-size triples in [1,4]^3 exclude contraction')
    print('SCOPE: exact finite checks only; the arbitrary-size family needs the general argument')


if __name__ == '__main__':
    fidelity()
    exact_witness()
    family()
    print('ALL_CHECKS_OK')
