#!/usr/bin/env python3
"""Exact checks for BGHS Question 3.21's boundary-dimension assertion.

Python 3 + SymPy. No floating point, random sampling, or enumeration of the
large witness is used. Finite checks support the stated general argument;
they are not a Lean proof or a verification of the separate degree question.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import comb, factorial
import json
import sympy as sp


def incidence(n, k, t):
    rows = tuple(map(frozenset, combinations(range(1, n + 1), t)))
    cols = tuple(map(frozenset, combinations(range(1, n + 1), k)))
    return rows, cols, sp.Matrix([[int(T <= F) for F in cols] for T in rows])


def image_exponents(monomial, t):
    return Counter(T for F in monomial for T in combinations(sorted(F), t))


def halves(pairs, base=()):
    result = [set(), set()]
    for bits in product((0, 1), repeat=len(pairs)):
        F = frozenset(base) | frozenset(pairs[i][bits[i]] for i in range(len(pairs)))
        assert len(F) == len(base) + len(pairs)
        result[sum(bits) % 2].add(F)
    assert not (result[0] & result[1])
    assert all(len(h) == 2 ** (len(pairs) - 1) for h in result)
    return result


def fidelity():
    # Definition 3.2's displayed A_(4,3,2), in exactly its printed row order.
    rows, cols, A = incidence(4, 3, 2)
    order = [(1, 2), (1, 3), (2, 3), (1, 4), (2, 4), (3, 4)]
    reordered = [list(A.row(rows.index(frozenset(T)))) for T in order]
    assert reordered == [[1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1],
                         [0, 1, 1, 0], [0, 1, 0, 1], [0, 0, 1, 1]]
    assert A.rank() == 4
    # Proposition 3.4 and Corollary 3.6 in two nonzero-kernel cases.
    for n, k, t in [(6, 3, 2), (8, 4, 3)]:
        _, cols, A = incidence(n, k, t)
        rank = A.rank()
        assert rank == min(comb(n, t), comb(n, k))
        assert len(cols) - rank == comb(n, k) - comb(n, t)
        print(json.dumps(dict(check='source-rank', n=n, k=k, t=t,
                              rank=rank, kernel_rank=len(cols)-rank,
                              projective_toric_dimension=rank-1), sort_keys=True))
    # Both explicit source Example 3.5 binomials lie in the monomial-map kernel.
    quartic = (['136', '246', '145', '235'], ['146', '236', '245', '135'])
    sextic = (['146', '156', '236', '123', '345', '245'],
              ['136', '126', '456', '145', '234', '235'])
    for degree, terms in [(4, quartic), (6, sextic)]:
        sides = [[frozenset(map(int, s)) for s in side] for side in terms]
        assert image_exponents(sides[0], 2) == image_exponents(sides[1], 2)
        print(json.dumps(dict(check='source-binomial', degree=degree, image_equal=True)))


def generator_checks():
    for k in range(2, 9):
        pairs = [(2*i, 2*i+1) for i in range(k)]
        E, O = halves(pairs)
        assert image_exponents(E, k-1) == image_exponents(O, k-1)
        swapped = [(pairs[0][1], pairs[0][0])] + pairs[1:]
        assert halves(swapped)[0] == O  # one tuple encoding includes both signs
        print(json.dumps(dict(check='parity-halves', k=k, facets=2**k,
                              half_cardinality=len(E), incidence_equal=True)))
    # Cone version: (t+1)-pairs with k-(t+1) common vertices.
    for t in range(1, 5):
        for cone_size in range(1, 4):
            pairs = [(2*i, 2*i+1) for i in range(t+1)]
            base = tuple(range(2*(t+1), 2*(t+1)+cone_size))
            E, O = halves(pairs, base)
            assert image_exponents(E, t) == image_exponents(O, t)
    print('GENERAL_POD_INCIDENCE_OK t=1..4 cone_size=1..3')


def finite_alteration_check():
    # Enumerate every sample for n=4,k=2, not the large proposed witness.
    W = tuple(map(frozenset, combinations(range(4), 2)))
    H = set()
    for tup in permutations(range(4)):
        E, O = halves([tup[:2], tup[2:]])
        H.add(frozenset(E))
        assert frozenset(O) in {frozenset(halves([u[:2], u[2:]])[0])
                                for u in permutations(range(4))}
    assert len(H) == 3 <= 4**4
    p = Q(1, 3)
    expectation = Q(0)
    best = 0
    for bits in product((0, 1), repeat=len(W)):
        chosen = {F for F, b in zip(W, bits) if b}
        bad = [e for e in H if e <= chosen]
        deleted = {min(e, key=lambda F: tuple(sorted(F))) for e in bad}
        survivor = chosen - deleted
        assert len(deleted) <= len(bad)
        assert not any(e <= survivor for e in H)
        assert len(survivor) >= len(chosen)-len(bad)
        expectation += p**len(chosen)*(1-p)**(len(W)-len(chosen))*(len(chosen)-len(bad))
        best = max(best, len(survivor))
    assert expectation == p*len(W)-len(H)*p**2
    assert best >= expectation
    print(json.dumps(dict(check='all-64-alteration-outcomes', edges=len(H),
                          expectation=str(expectation), best_survivor=best), sort_keys=True))


def counterfamily_checks():
    assert Q(1, 192)-Q(1, 65536) > Q(1, 200)
    for a in [4, 5, 8, 16]:
        n, p = a**7, Q(1, 4*a**4)
        M, B = comb(n, 4), n**8
        lower = p*M-B*p**8
        threshold = comb(n, 3)
        assert 0 < p <= 1
        assert M >= Q(n**4, 48)
        assert lower >= a**24*(Q(1, 192)-Q(1, 65536))
        assert a**24*(Q(1, 192)-Q(1, 65536)) > Q(a**24, 200)
        assert Q(a**24, 200) > Q(a**21, 6) >= threshold
        assert lower > threshold
        guaranteed_cardinality = -(-lower.numerator // lower.denominator)
        print(json.dumps(dict(check='counterfamily-bound', a=a, n=n, p=str(p),
                              variables=M, edge_count_upper=B,
                              alteration_lower=str(lower),
                              guaranteed_coordinate_dimension=guaranteed_cardinality-1,
                              toric_dimension=threshold-1,
                              strict_dimension_gap=True), sort_keys=True))
    # Check factored constants/exponents in the general k argument.
    for k in range(4, 13):
        r, C = 2**(k-1), 2**(k+1)*factorial(k)
        a = C*C
        assert r-k-1 > 0
        assert Q(2, C*C)-Q(1, C**r) >= Q(1, C*C)
        assert k*(r-2)-(r-1)*(k-1) == r-k-1
        assert a**(r-k-1)*factorial(k-1) > C*C
        print(json.dumps(dict(check='general-family-constants', k=k, r=r, C=C,
                              exponent_gap=r-k-1), sort_keys=True))


if __name__ == '__main__':
    fidelity()
    generator_checks()
    finite_alteration_check()
    counterfamily_checks()
    print('ALL_CHECKS_OK')
