#!/usr/bin/env python3
"""Exact checks for Ben Shimon--Zabokritskiy, arXiv:2607.19566v1, C15.

Independent implementation from the paper's definitions. Standard library only;
no companion-package program source was read or used. Positions in the printed
witnesses are one-based. Assertions are required: do not run with python -O.
"""
from functools import lru_cache
from itertools import permutations, product
from math import gcd, lcm


@lru_cache(maxsize=None)
def swap_sets(n):
    """All matchings of the path on n positions, as zero-based swap sites."""
    if n <= 1:
        return ((),)
    without_last = swap_sets(n - 1)
    with_last = tuple(s + (n - 2,) for s in swap_sets(n - 2))
    return without_last + with_last


def apply_swaps(word, sites):
    assert all(0 <= i < len(word) - 1 for i in sites)
    assert all(j - i >= 2 for i, j in zip(sites, sites[1:]))
    out = list(word)
    for i in sites:
        out[i], out[i + 1] = out[i + 1], out[i]
    return ''.join(out)


@lru_cache(maxsize=None)
def truncated_outputs(word):
    return frozenset(apply_swaps(word, s)[:-1] for s in swap_sets(len(word)))


def prefixes(first, blocks, n):
    """Enumerate whole-block extensions up to their first crossing of n.

    Every longer concatenation has the same n-prefix as its shortest block
    prefix reaching n; nonempty blocks make this finite search exhaustive.
    """
    assert n >= 1 and all(blocks)
    frontier = {first}
    result = set()
    while frontier:
        word = frontier.pop()
        if len(word) >= n:
            result.add(word[:n])
        else:
            frontier.update(word + block for block in blocks)
    return frozenset(result)


def output_union(first, blocks, n):
    return frozenset().union(*(truncated_outputs(w) for w in prefixes(first, blocks, n)))


def common_outputs(a, b, n):
    blocks = (a, b)
    return output_union(a, blocks, n) & output_union(b, blocks, n)


def first_extinction(a, b, limit):
    for n in range(1, limit + 1):
        if not common_outputs(a, b, n):
            return n
    return None


def permutation_ball(word):
    """Separate literal radius-one displacement-permutation implementation."""
    return frozenset(''.join(word[p[i]] for i in range(len(word)))
                     for p in permutations(range(len(word)))
                     if all(abs(p[i] - i) <= 1 for i in range(len(word))))


def check_channel_fidelity():
    cases = 0
    for n in range(1, 7):
        for letters in product('01', repeat=n):
            w = ''.join(letters)
            swaps = frozenset(apply_swaps(w, s) for s in swap_sets(n))
            assert swaps == permutation_ball(w), (n, w)
            cases += 1
    assert len(swap_sets(11)) == 144
    print(f'Channel fidelity: {cases} binary words of lengths 1..6; displacement ball = matching ball')
    print('Legal swap sets at length 11: 144')


def check_source_results():
    assert first_extinction('00', '11', 2) == 2
    print('Source first-extinction example: K*(00,11) = 2')
    for length in range(3, 9):
        a, b = '0' * length, '0' * (length - 1) + '1'
        depth = first_extinction(a, b, length + 2)
        assert depth == length + 2, (length, depth)
        assert first_extinction(b, a, length + 2) == depth
        print(f'Source Proposition 12 fidelity: L={length}, K*={depth}=L+2 (both orientations)')


def one_positions(word):
    return tuple(i + 1 for i, c in enumerate(word) if c == '1')


def check_counterexample():
    a, b = '001', '100000'
    assert a and b and a != b and not a.startswith(b) and not b.startswith(a)
    assert (len(a), len(b), gcd(len(a), len(b)), lcm(len(a), len(b))) == (3, 6, 3, 6)
    ceiling = lcm(len(a), len(b)) + len(a) - 1
    assert ceiling == 8
    counts = [len(common_outputs(a, b, n)) for n in range(1, 12)]
    assert counts == [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 0], counts
    print('Candidate a=001, b=100000: ordered lengths=3,6; gcd=3; lcm=6; proposed ceiling=8')
    print('Common-output counts at depths 1..11:', ','.join(map(str, counts)))
    assert first_extinction(a, b, 11) == first_extinction(b, a, 11) == 11
    u, v = '0010010010', '1000001000'
    assert u == (a * 4)[:10] and v == (b * 2)[:10]
    assert u in prefixes(a, (a, b), 10) and v in prefixes(b, (a, b), 10)
    left = apply_swaps(u, (1, 5, 8))
    right = apply_swaps(v, (0,))
    assert (left, right) == ('0100001001', '0100001000')
    assert left[:-1] == right[:-1] == '010000100'
    print('Depth-10 witness:', u, 'sites={2,6,9} ->', left, ';', v, 'sites={1} ->', right)
    for m in range(1, 11):
        ls = tuple(i for i in (1, 5, 8) if i + 1 < m)
        rs = tuple(i for i in (0,) if i + 1 < m)
        assert apply_swaps(u[:m], ls)[:-1] == apply_swaps(v[:m], rs)[:-1]
    left_inputs = prefixes(a, (a, b), 11)
    right_inputs = prefixes(b, (a, b), 11)
    assert left_inputs == frozenset(['00100100100', '00100100110', '00100110000',
                                    '00110000000', '00110000010'])
    assert right_inputs == frozenset(['10000000100', '10000000110', '10000010000'])
    print('Complete left depth-11 inputs:', ','.join(sorted(left_inputs)))
    print('Complete right depth-11 inputs:', ','.join(sorted(right_inputs)))
    left_outputs = output_union(a, (a, b), 11)
    right_outputs = output_union(b, (a, b), 11)
    right_zero = {w for w in right_outputs if w[0] == '0'}
    assert right_zero == {'0100000001', '0100000010', '0100000011', '0100000100',
                          '0100000101', '0100001000', '0100010000'}
    print('Right truncated outputs beginning 0:', ','.join(sorted(right_zero)))
    assert all(w[0] == '0' for w in left_outputs)
    for w in right_zero:
        positions = one_positions(w)
        assert (len(positions) == 2 and positions[1] >= 6) or (
                len(positions) == 3 and positions[1] >= 8)
    for w in sorted(left_inputs):
        outputs01 = {z for z in truncated_outputs(w) if z.startswith('01')}
        assert outputs01
        for s in swap_sets(11):
            z = apply_swaps(w, s)[:-1]
            if z.startswith('01'):
                assert 1 in s and 0 not in s and 2 not in s
        if w in {'00100100100', '00100100110', '00100110000'}:
            assert all(len(one_positions(z)) >= 3 and one_positions(z)[1] <= 7 for z in outputs01)
        else:
            assert all(len(one_positions(z)) in {2, 3} and one_positions(z)[1] in {4, 5} for z in outputs01)
        signatures = sorted({(len(one_positions(z)), one_positions(z)[1]) for z in outputs01})
        print('Left input', w, 'starting-01 output (one-count, second-one-position):', signatures)
    assert all(not (truncated_outputs(u11) & truncated_outputs(v11))
               for u11 in left_inputs for v11 in right_inputs)
    assert not (left_outputs & right_outputs)
    print('All 5 x 3 depth-11 prefix pairs extinct; positional separation verified on every output')
    print('Expected refutation confirmed by exact enumeration: finite K*=11 > 8')


def main():
    if not __debug__:
        raise RuntimeError('Assertions must be enabled')
    check_channel_fidelity()
    check_source_results()
    check_counterexample()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
