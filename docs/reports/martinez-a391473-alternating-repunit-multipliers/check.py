#!/usr/bin/env python3
"""Exact, independent finite checks for both OEIS A391473 conjectures.

No external packages, network, random choices, or floating-point arithmetic.
The unbounded minimality assertion is not established by these checks.
"""
from collections import deque
from fractions import Fraction


OEIS_TERMS = (
    1, 11, 19, 371, 1891, 36911, 189091, 3690911, 18909091,
    369090911, 1890909091, 36909090911, 189090909091,
    3690909090911, 18909090909091, 369090909090911,
    1890909090909091,
)


def repunit(length):
    assert length >= 0
    numerator = 10**length - 1
    assert numerator % 9 == 0
    return numerator // 9


def alternating(value):
    """Canonical positive decimal expansion, examined from units upward."""
    if value <= 0:
        return False
    previous = value % 2
    value //= 10
    while value:
        parity = (value % 10) % 2
        if parity == previous:
            return False
        previous = parity
        value //= 10
    return True


def formula(length):
    assert length >= 3
    if length % 2:
        numerator = 208 * 10**(length - 3) + 1
    else:
        numerator = 4060 * 10**(length - 4) + 21
    quotient, remainder = divmod(numerator, 11)
    assert remainder == 0
    assert quotient > 0
    return quotient


def brute_force(length):
    """Enumerate multipliers in increasing order; do not use the conjecture."""
    r = repunit(length)
    multiplier = 1
    while not alternating(multiplier * r):
        multiplier += 1
    return multiplier


def automaton_minimum(length):
    """Breadth-first digit search, retaining the first prefix for each state.

    A state is (remainder modulo the repunit, last-digit parity). Equal
    states have identical allowed suffixes. Ordered BFS finds the shortest,
    then lexicographically first positive alternating divisible integer.
    """
    modulus = repunit(length)
    parents = {}
    queue = deque()

    def reconstruct(state):
        digits = []
        while state is not None:
            state, digit = parents[state]
            digits.append(str(digit))
        return int(''.join(reversed(digits)))

    for digit in range(1, 10):
        state = (digit % modulus, digit % 2)
        if state in parents:
            continue
        parents[state] = (None, digit)
        if state[0] == 0:
            return reconstruct(state)
        queue.append(state)
    while queue:
        state = queue.popleft()
        residue, last_parity = state
        for digit in range(1 - last_parity, 10, 2):
            following = ((10 * residue + digit) % modulus, digit % 2)
            if following in parents:
                continue
            parents[following] = (state, digit)
            if following[0] == 0:
                return reconstruct(following)
            queue.append(following)
    raise AssertionError('no alternating multiple found')


def check_blocks(length):
    h = length - 2
    t = 10**h
    if h % 2:
        q = (h - 1) // 2
        high = 2 * 100**q + 10 * sum(100**i for i in range(q))
        j = 1
        expected_high = '2' + '10' * q
        expected_low = '9' + '01' * q
        high_closed = Fraction(208 * 10**(h - 1) - 10, 99)
    else:
        q = h // 2
        high = 4 * 10**(h - 1) + sum(100**i for i in range(q))
        j = 2
        expected_high = '4' + '10' * (q - 1) + '1'
        expected_low = '8' + '12' * (q - 1) + '1'
        high_closed = Fraction(406 * 10**(h - 1) - 1, 99)
    low = t + j * repunit(h) - high
    assert high == high_closed
    assert str(high) == expected_high
    assert str(low).zfill(h) == expected_low
    multiplier = 9 * high + j
    multiple = multiplier * repunit(length)
    assert multiplier == formula(length)
    assert str(multiple) == expected_high + str(j) + str(j - 1) + expected_low
    assert alternating(multiple)
    assert high + low == t + j * repunit(h)
    assert 10**(2 * length - 3) <= multiple < 10**(2 * length - 2)
    high_check, tail = divmod(multiple, 10**length)
    assert high_check == high
    assert tail == j * repunit(length) - high
    if length % 2 == 0:
        high_digits = [int(d) for d in reversed(str(high))]
        low_digits = [int(d) for d in reversed(str(low).zfill(h))]
        carry = 0
        for index, (d, e) in enumerate(zip(high_digits, low_digits)):
            assert d % 2 == e % 2
            assert carry == 0
            carry, residue = divmod(d + e + carry, 10)
            assert residue == j
            assert carry == (1 if index == h - 1 else 0)
        assert carry == 1


def main():
    assert not alternating(0)
    assert alternating(1) and alternating(121) and alternating(1234567890)
    assert not alternating(11) and not alternating(20)
    print('Fidelity: brute-force least positive multipliers (no formula used)')
    print('length repunit multiplier product')
    for length in range(1, 9):
        multiplier = brute_force(length)
        assert multiplier == OEIS_TERMS[length - 1]
        print(length, repunit(length), multiplier, multiplier * repunit(length))
    print('Independent remainder/parity BFS, lengths 1..6:')
    for length in range(1, 7):
        multiple = automaton_minimum(length)
        assert multiple % repunit(length) == 0
        assert alternating(multiple)
        multiplier = multiple // repunit(length)
        assert multiplier == OEIS_TERMS[length - 1]
        print(length, multiplier, multiple)
    for length in range(3, len(OEIS_TERMS) + 1):
        assert formula(length) == OEIS_TERMS[length - 1]
    print('Both formulas agree with every published term of lengths 3..17.')
    for length in range(3, 201):
        check_blocks(length)
    print('Exact construction, geometric sums, digit patterns, and carry checks: lengths 3..200.')
    print('Finite checks support the proposed proof; universal minimality requires the carry argument.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
