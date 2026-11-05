#!/usr/bin/env python3
"""Exact checks for Alm--Maddux arXiv:1712.00129v1, Problem 3.

No external packages. Integers encode bitstrings, bit i is coordinate i;
the last coordinate is bit 3*k. All convolutions use integer arithmetic.
The source's 63 bitstrings below are transcribed literally, parsed MSB first.
"""
from fractions import Fraction
from itertools import permutations
from math import comb
import json

SOURCE_B_WORDS = """0110000000 1000000010 0000000101 1000000111 0100001000 1100001010 0100001101 1100001111 0000010001 1000010011 0000010100 1000010110 0010011001 0100011001 1010011011 0100011100 1100011110 0110000101 1000100000 0000100010 1100011011 1000100101 0000100111 1100101000 0010011100 0100101010 1110000111 1100101101 0100101111 1000110001 0000110011 1000110100 1010011110 0000110110 1100111001 0100111011 1100111100 1010001010 0100111110 1110100000 0010101111 0110110110 1110110100 0110110011 1110010110 0110100010 0010001101 1010111001 1010001111 1110000010 1110100101 1010101101 0010111011 1110110001 0110010001 1010111100 0110100111 1010101000 1110010011 0010111110 0110010100 0010001000 0010101010"""


def direct_sum(left, right):
    return {x ^ y for x in left for y in right}


def fwht(values):
    out = list(values)
    half = 1
    while half < len(out):
        for start in range(0, len(out), 2 * half):
            for i in range(start, start + half):
                x, y = out[i], out[i + half]
                out[i], out[i + half] = x + y, x - y
        half *= 2
    return out


def convolution(left_transform, right_transform):
    raw = fwht([x * y for x, y in zip(left_transform, right_transform)])
    size = len(raw)
    assert all(x % size == 0 for x in raw)
    result = [x // size for x in raw]
    assert all(x >= 0 for x in result)
    return result


def required_sets(a, b, c, group):
    return {'AA': group, 'AB': a | c, 'AC': group - {0},
            'BB': b | {0}, 'BC': a | c, 'CC': group - c}


def partition_check(a, b, c, size):
    assert a and b and c
    assert not (a & b or a & c or b & c)
    assert a | b | c == set(range(1, size))


def table_from_source_cycles():
    # Independent cycle-to-table check, including all permutations and identities.
    atoms = 'iabc'
    allowed = set()
    for cycle in ['iii', 'iaa', 'ibb', 'icc',
                  'aaa', 'bbb', 'acc', 'aab', 'aac', 'bcc', 'abc']:
        allowed.update(permutations(cycle))
    table = {(i, j): {z for z in atoms if (i, j, z) in allowed}
             for i in atoms for j in atoms}
    expected = {'aa': set('iabc'), 'ab': set('ac'), 'ac': set('abc'),
                'bb': set('ib'), 'bc': set('ac'), 'cc': set('iab')}
    for pair, products in expected.items():
        assert table[tuple(pair)] == products
        assert table[tuple(reversed(pair))] == products
    for i in atoms:
        assert table['i', i] == table[i, 'i'] == {i}
    for cycle in ['ccc', 'abb', 'cbb']:
        assert all(tuple(t) not in allowed for t in permutations(cycle))
    return {'mandatory_cycles': ['aaa', 'bbb', 'acc', 'aab', 'aac', 'bcc', 'abc'],
            'forbidden_cycles': ['ccc', 'abb', 'cbb'], 'table_check': 'OK'}


def source_fidelity():
    # Source Section 2: its explicit dimension-ten subgroup, not the new construction.
    size = 1 << 10
    group = set(range(size))
    b = {int(word, 2) for word in SOURCE_B_WORDS.split()}
    h = b | {0}
    x = {g for g in group if 1 <= g.bit_count() <= 6}
    c = {g for g in group if g.bit_count() > 6}
    a = x - b
    assert len(b) == 63 and len(h) == 64
    assert direct_sum(h, h) == h and h <= x | {0}
    partition_check(a, b, c, size)
    assert direct_sum(x, x) == group
    assert direct_sum(x, c) == group - {0}
    assert direct_sum(c, c) == group - c
    sets = {'A': a, 'B': b, 'C': c}
    for pair, expected in required_sets(a, b, c, group).items():
        assert direct_sum(sets[pair[0]], sets[pair[1]]) == expected
    n = 13
    points = comb(3 * n - 4, n)
    bound = points ** 2 * 4 ** 3 * Fraction(2, 3) ** ((n - 2) ** 2)
    assert points == 1476337800 and bound < 1
    return {'dimension': 10, 'source_bitstrings': len(b),
            'atom_sizes': [len(a), len(b), len(c)],
            'six_direct_sumsets': 'OK', 'johnson_n': n, 'johnson_points': points,
            'johnson_union_bound_less_than_one': True}


def family_check(k):
    dimension = 3 * k + 1
    size = 1 << dimension
    group = set(range(size))
    h = {g for g in group if g >> (3 * k) == 0 and
         all(((g >> (3 * j)) & 7).bit_count() % 2 == 0 for j in range(k))}
    b = h - {0}
    c = {g for g in group if g.bit_count() > 2 * k}
    a = group - h - c
    partition_check(a, b, c, size)
    assert len(h) == 4 ** k and all(g.bit_count() <= 2 * k for g in h)
    assert direct_sum(h, h) == h
    sets = {'A': a, 'B': b, 'C': c}
    transforms = {name: fwht([int(g in part) for g in range(size)])
                  for name, part in sets.items()}
    readings = {}
    for pair, expected in required_sets(a, b, c, group).items():
        counts = convolution(transforms[pair[0]], transforms[pair[1]])
        actual = {g for g, count in enumerate(counts) if count > 0}
        missing, excess = sorted(expected - actual), sorted(actual - expected)
        if k == 2 and pair == 'BC':
            assert missing == [63] and excess == []
        else:
            assert missing == [] and excess == []
        if k == 3:
            assert actual == direct_sum(sets[pair[0]], sets[pair[1]])
        readings[pair] = {'missing': missing, 'excess': excess,
                          'min_required_witness_count': min(counts[g] for g in expected)}
    return {'k': k, 'dimension': dimension, 'group_size': size,
            'atom_sizes': [len(a), len(b), len(c)], 'sumsets': readings,
            'direct_convolution_crosscheck': k == 3}


def binomial_check():
    # Finite checks of the stated recurrence and strict bound, not a universal proof.
    for k in range(2, 101):
        n = 3 * k + 1
        s = sum(comb(n, i) for i in range(k + 1))
        next_s = sum(comb(n + 3, i) for i in range(k + 2))
        assert next_s == 8 * s + comb(n, k + 1) - 4 * comb(n, k) - comb(n, k - 1)
        assert Fraction(comb(n, k + 1), comb(n, k)) == Fraction(2 * k + 1, k + 1)
        assert next_s < 8 * s and s + 4 ** k < 8 ** k
    return {'integer_parameters': [2, 100], 'recurrence_and_strict_bound': 'OK'}


def main():
    print(json.dumps({'source_cycle_table': table_from_source_cycles()}, sort_keys=True))
    print(json.dumps({'source_fidelity': source_fidelity()}, sort_keys=True))
    for k in range(2, 6):
        print(json.dumps({'family': family_check(k)}, sort_keys=True))
    print(json.dumps({'binomial': binomial_check()}, sort_keys=True))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
