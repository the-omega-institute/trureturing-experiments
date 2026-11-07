#!/usr/bin/env python3
"""Exact finite checks for Zhu--Richard, arXiv:2610.10967v1, Section 4.2.

The full functional graph is enumerated, without the proposed reduction, for
all source-defined sizes 4 <= n <= 22. Kahn deletion leaves precisely cycles.
Only integers and Fraction are used as evidence. This is not an unbounded proof.
"""
from array import array
from fractions import Fraction


def local(a, b, c):
    return c ^ (a & (1 ^ b))


def step(x, n, rule=154):
    mask = (1 << n) - 1
    left = ((x << 1) & mask) | (x >> (n - 1))
    right = (x >> 1) | ((x & 1) << (n - 1))
    if rule == 154:
        return right ^ (left & (mask ^ x))
    return sum(((rule >> (4 * ((left >> i) & 1) +
                            2 * ((x >> i) & 1) + ((right >> i) & 1))) & 1) << i
               for i in range(n))


def cyclic_count(successors):
    indegree = array('I', [0]) * len(successors)
    for y in successors:
        indegree[y] += 1
    queue = array('I', (x for x, degree in enumerate(indegree) if degree == 0))
    head = 0
    while head < len(queue):
        x = queue[head]
        head += 1
        y = successors[x]
        indegree[y] -= 1
        if indegree[y] == 0:
            queue.append(y)
    return len(successors) - len(queue)


def transient_formula(n):
    if n % 2:
        return 0
    m = n // 2
    d = m & -m
    return 2 * ((1 << m) - (1 << (m - d)))


def parity_sector(x, n):
    even = sum(1 << i for i in range(0, n, 2))
    odd = ((1 << n) - 1) ^ even
    return x & even == 0 or x & odd == 0


def linear_step(a, m):
    mask = (1 << m) - 1
    return (((a << 1) & mask) | (a >> (m - 1))) ^ ((a >> 1) | ((a & 1) << (m - 1)))


def main():
    truth_table = [local(a, b, c) for a in range(2) for b in range(2) for c in range(2)]
    assert truth_table == [0, 1, 0, 1, 1, 0, 0, 1]
    assert sum(bit << i for i, bit in enumerate(truth_table)) == 154
    print('WOLFRAM_TRUTH_TABLE', truth_table)
    # Source Example 2.1/Figure 3: Rule 74 at n=4 has one fixed point and a 4-cycle.
    s74 = [step(x, 4, 74) for x in range(16)]
    assert cyclic_count(s74) == 5
    assert sum(x == y for x, y in enumerate(s74)) == 1
    print('SOURCE_FIDELITY Rule74 n=4 states=16 periodic=5 nonperiodic=11 fixed=1')
    ratios = {}
    for n in range(4, 23):
        succ = array('I', (step(x, n) for x in range(1 << n)))
        periodic = cyclic_count(succ)
        transient = (1 << n) - periodic
        expected = transient_formula(n)
        assert transient == expected, (n, transient, expected)
        z = Fraction(transient, 1 << n)
        ratios[n] = z
        if n >= 13:
            assert z < Fraction(1, 100)
        print(f'FULL_GRAPH n={n} states={1 << n} periodic={periodic} nonperiodic={transient} z={z}')
        if n <= 12:
            fibers = {}
            for x, y in enumerate(succ):
                fibers.setdefault(y, []).append(x)
            for fiber in fibers.values():
                if len(fiber) > 1:
                    assert n % 2 == 0
                    assert all(parity_sector(x, n) for x in fiber), (n, fiber)
            if n % 2 == 0:
                assert all(parity_sector(x, n) == parity_sector(succ[x], n)
                           for x in range(1 << n))
    print('COLLISION_AND_PARITY_INVARIANCE_CHECK sizes=4..12 OK')
    # Independently enumerate the proposed compressed linear graph, not its formula.
    for m in range(2, 17):
        lin = [linear_step(a, m) for a in range(1 << m)]
        periodic = cyclic_count(lin)
        d = m & -m
        assert periodic == 1 << (m - d), (m, periodic)
        n = 2 * m
        sector_transient = 2 * ((1 << m) - periodic)
        assert sector_transient == transient_formula(n)
        for a in range(1 << m):
            x = sum(((a >> j) & 1) << (2 * j) for j in range(m))
            y = sum(((lin[a] >> j) & 1) << (2 * j) for j in range(m))
            assert step(step(x, n), n) == y, (m, a)
        print(f'COMPRESSED_GRAPH n={n} m={m} periodic_per_sector={periodic} nonperiodic_in_union={sector_transient}')
    for n in range(13, 28):
        z = Fraction(transient_formula(n), 1 << n)
        assert z < Fraction(1, 100)
        print(f'SOURCE_RANGE n={n} formula_z={z} below_1/100=True')
    # Figure 13 displays rounded-to-eight-decimal annotations, checked exactly.
    for n, printed in [(22, Fraction('0.00048828')), (24, Fraction('0.00045776'))]:
        z = Fraction(transient_formula(n), 1 << n)
        assert abs(z - printed) <= Fraction(1, 2 * 10**8)
        print(f'SOURCE_FIGURE13 n={n} z={z} printed_decimal_as_fraction={printed} rounding_error={abs(z-printed)}')
    previous = Fraction(transient_formula(4), 1 << 4)
    for m in range(3, 257):
        z = Fraction(transient_formula(2 * m), 1 << (2 * m))
        assert z < previous
        assert Fraction(1, 1 << m) <= z < Fraction(1, 1 << (m - 1))
        previous = z
    print('EXACT_MONOTONICITY_AND_ENVELOPE_CHECK even_sizes=4..512 OK')
    for n, expected in [(28, 24576), (30, 32768), (32, 131070)]:
        assert transient_formula(n) == expected
        print(f'BEYOND_SOURCE n={n} nonperiodic={expected} z={Fraction(expected,1 << n)}')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
