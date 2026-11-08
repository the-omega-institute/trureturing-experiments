#!/usr/bin/env python3
"""Exact checks for HKLT arXiv:2605.18644v1, Section 2, Problem 4.

Standard library only. This program supplies experimental evidence, not a
Lean proof. Congruence pairs are (modulus, residue).
"""
from fractions import Fraction as Q
from itertools import combinations, product
from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sums(moduli):
    s1 = sum((Q(1, m) for m in moduli), Q())
    s2 = sum((Q(1, m*n) for m, n in combinations(moduli, 2)
              if gcd(m, n) == 1), Q())
    s3 = sum((Q(1, m*n*k) for m, n, k in combinations(moduli, 3)
              if gcd(m, n) == gcd(n, k) == gcd(k, m) == 1), Q())
    return s1, s2, s3


def covered(system, x):
    return any((x-r) % m == 0 for m, r in system)


def density(system):
    period = lcm(*(m for m, _ in system))
    return Q(sum(covered(system, x) for x in range(1, period+1)), period)


def block(cutoff):
    return [m for m in range(3, cutoff+1, 6) if m % 11]


def source_fidelity():
    # Introduction, equation (1.1): Davenport's displayed cover, m=2, M=L=12.
    davenport = [(2, 1), (4, 2), (3, 0), (6, 4), (12, 8)]
    require(lcm(*(m for m, _ in davenport)) == 12, 'source period')
    require(density(davenport) == 1, 'source Davenport cover')
    s1, s2, s3 = sums([m for m, _ in davenport])
    require((s1, s2, s3) == (Q(4, 3), Q(1, 4), Q()), 'source sums')
    require(density(davenport) <= s1-s2+s3, 'source Theorem 1.9 example')
    # The sentence after Problem 4 says the LCM condition can be dropped at n=3.
    cases = 0
    for moduli in combinations(range(2, 9), 3):
        bound = sum((a*b for a, b in zip(sums(moduli), [1, -1, 1])), Q())
        period = lcm(*moduli)
        # Integer masks count [1,L] without rational approximations.
        masks = [[sum(1 << (x-1) for x in range(1, period+1) if x % m == r)
                  for r in range(m)] for m in moduli]
        for rs in product(*(range(m) for m in moduli)):
            mask = masks[0][rs[0]] | masks[1][rs[1]] | masks[2][rs[2]]
            require(Q(mask.bit_count(), period) <= bound, 'source n=3 assertion')
            cases += 1
    print('SOURCE_DAVENPORT density=1 period=12 bound=13/12')
    print('SOURCE_N3 moduli=2..8 all_residue_assignments=' + str(cases))


def counterexamples():
    a_system = [(2, 0), (4, 1), (8, 3), (16, 7), (32, 15), (64, 31),
                (128, 63), (256, 127), (512, 255), (11, 0), (22, 1),
                (44, 35), (88, 47), (176, 15), (352, 159), (704, 127),
                (1408, 1151), (2816, 767), (5632, 4607)]
    require(a_system[:9] == [(2**j, 2**(j-1)-1) for j in range(1, 10)],
            'binary classes use exact valuation of x+1')
    for j, (m, r) in enumerate(a_system[9:]):
        require(m == 11*2**j and (r+1) % 2**j == 0 and r % 11 == j,
                'CRT residue table')
    require([x for x in range(5632) if not covered(a_system, x)] == [5631],
            'A misses exactly -1 mod 5632')
    a, e, t = sums([m for m, _ in a_system])
    require((a, e, t) == (Q(151, 128), Q(511, 5632), Q()), 'A exact sums')
    require(1-a+e == -Q(501, 5632), 'negative coefficient')
    b489 = sum((Q(1, m) for m in block(489)), Q())
    require(len(block(489)) == 75, '75 added moduli')
    require(b489 == Q(3724374094737398357076456355353489210538987739102697298950471495779,
                      3715236714461777943939088897767713808889287200964363185634535244625),
            'exact reciprocal certificate')
    floors = sum(1000000 // m for m in block(489))
    require(floors == 1002429 and b489 >= Q(floors, 1000000) > Q(501, 500),
            'independent floor certificate')
    reduced_period = 16896
    for cutoff in [489, 600, 1000]:
        bs = block(cutoff)
        system = a_system + [(m, 0) for m in bs]
        ms = [m for m, _ in system]
        require(len(ms) == len(set(ms)) and min(ms) > 1, 'distinct nonunit moduli')
        require(3 in bs and all(m % 3 == 0 and gcd(m, 22) == 1 for m in bs),
                'B union is exactly 0 mod 3, all cross pairs coprime')
        require(all(gcd(m, n) > 1 for m, n in combinations(bs, 2)), 'B has no edges')
        b = sum((Q(1, m) for m in bs), Q())
        s1, s2, s3 = sums(ms)
        require((s1, s2, s3) == (a+b, e+a*b, e*b), 'independent pair/triple enumeration')
        upper = s1-s2+s3
        require(upper == 1-Q(501, 5632)*(b-1), 'family identity')
        missed = [x for x in range(reduced_period) if not covered(system, x)]
        require(missed == [11263, 16895], 'reduced-period uncovered residues')
        actual = Q(reduced_period-len(missed), reduced_period)
        require(actual == Q(8447, 8448) and upper < Q(5631, 5632) < actual,
                'strict counterexample')
        require(not covered(system, -1), 'negative integer uncovered')
        print(f'FAMILY cutoff={cutoff} moduli={len(ms)} density={actual} bound={upper}')
        if cutoff == 489:
            print('WITNESS_MAX_MODULUS', max(ms))
            print('WITNESS_LCM', lcm(*ms))
            print('WITNESS_B_MASS', b)
            print('FLOOR_CERTIFICATE', floors)
            print('UNCOVERED_MOD_16896', missed)
    residues = [1, 3, 5, 7, 9, 13, 15, 17, 19, 21]
    for t in range(1000):
        ms = [3*(22*t+r) for r in residues]
        require(len(set(ms)) == 10 and all(m >= 3 and m % 6 == 3 and m % 11 for m in ms),
                'harmonic-block eligibility')
        require(sum((Q(1, m) for m in ms), Q()) >= Q(10, 3*(22*t+21)) >= Q(5, 33*(t+1)),
                'harmonic lower bound')
    print('DIVERGENCE_BLOCK_CHECK t=0..999 coefficient=-501/5632')


def chordal_pointwise_check():
    cases = 0
    for n in range(6):
        possible = list(combinations(range(n), 2))
        for bits in range(1 << len(possible)):
            edges = {pair for i, pair in enumerate(possible) if bits >> i & 1}
            if any(any((j, k) not in edges for j, k in combinations(
                    [j for j in range(i+1, n) if (i, j) in edges], 2)) for i in range(n)):
                continue
            for subset in range(1 << n):
                vs = [i for i in range(n) if subset >> i & 1]
                ne = sum((i, j) in edges for i, j in combinations(vs, 2))
                nt = sum(all(e in edges for e in combinations(triple, 2))
                         for triple in combinations(vs, 3))
                require(len(vs)-ne+nt >= int(bool(vs)), 'PEO pointwise inequality')
                cases += 1
    print('PEO_GRAPH_POINTWISE_CHECK n=0..5 cases=' + str(cases))


if __name__ == '__main__':
    source_fidelity()
    counterexamples()
    chordal_pointwise_check()
    print('ALL_CHECKS_OK')
