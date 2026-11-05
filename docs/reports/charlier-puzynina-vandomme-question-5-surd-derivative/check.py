#!/usr/bin/env python3
"""Exact finite controls for CPV arXiv:1907.00192v2, Question 5.

No dependencies beyond Python 3's standard library. Finite controls are not a
proof of the universally quantified counterexample family.
"""
from functools import lru_cache
from math import factorial, gcd, prod


def encode_first_appearance(words, start=0):
    table = {}
    codes = []
    for word in words:
        if word not in table:
            table[word] = len(table) + start
        codes.append(table[word])
    return tuple(codes), tuple(table)


def source_fidelity():
    # Section 8: return words to 011 in the Thue--Morse word and its
    # derivative prefix. Here the source uses one-based labels.
    bits = tuple(n.bit_count() % 2 for n in range(1 << 15))
    fixed = (0,)
    for _ in range(15):
        fixed = tuple(y for x in fixed for y in (x, 1 - x))
    assert fixed == bits
    starts = [i for i in range(len(bits) - 2) if bits[i:i + 3] == (0, 1, 1)]
    words = tuple(bits[a:b] for a, b in zip(starts, starts[1:]))
    codes, table = encode_first_appearance(words, start=1)
    expected_words = {'011010', '011001', '01101001', '0110'}
    assert {''.join(map(str, word)) for word in table} == expected_words
    expected_prefix = '12341243123431241234124'
    assert ''.join(map(str, codes[:len(expected_prefix)])) == expected_prefix
    print('FIDELITY: Thue-Morse return-word set = ' + ', '.join(sorted(expected_words)))
    print('FIDELITY: derivative prefix = ' + expected_prefix)


def word(m, x, y):
    return int((x + y) % m == m - 1)


def block(m, x, y, shape):
    return tuple(word(m, x + u, y + v)
                 for u in range(shape[0]) for v in range(shape[1]))


@lru_cache(maxsize=None)
def return_cycle(m, qx, qy, shape=(1, 1)):
    # Construct the sampled word literally for a full lattice period m.
    # The extra endpoint is its first letter again. Consecutive starts give
    # every return word in a cycle; repeat the cycle to check the boundary.
    samples = tuple(block(m, n * qx, n * qy, shape) for n in range(2 * m + 1))
    starts = [n for n, value in enumerate(samples) if value == samples[0]]
    words = tuple(samples[a:b] for a, b in zip(starts, starts[1:]))
    period_starts = sum(1 for n in starts if n < m)
    assert words[:period_starts] == words[period_starts:2 * period_starts]
    codes, table = encode_first_appearance(words)
    assert codes[:period_starts] == codes[period_starts:2 * period_starts]
    return codes[:period_starts], table


def derivative(m, x, y):
    if x == 0 and y == 0:
        return 0  # First-appearance code at index zero on every ray.
    radial = gcd(x, y)
    qx, qy = x // radial, y // radial
    assert gcd(qx, qy) == 1
    codes, _ = return_cycle(m, qx % m, qy % m)
    return codes[radial % len(codes)]


def finite_periodicity_and_rays():
    n_blocks = n_rays = 0
    for m in range(4, 13):
        # Exhaust all origins and direction residues, plus rectangular
        # shapes. Translation by m*q preserves the entire block.
        for x in range(m):
            for y in range(m):
                for qx in range(m):
                    for qy in range(m):
                        for shape in ((1, 1), (2, 1), (1, 2), (2, 3)):
                            assert block(m, x, y, shape) == block(
                                m, x + m * qx, y + m * qy, shape)
                            n_blocks += 1
        # Check the uniform-window arithmetic without enumerating origins.
        for k in range(3 * m):
            n = ((k + m - 1) // m) * m
            assert k <= n < k + m
        row = tuple(word(m, x, 0) for x in range(m))
        col = tuple(word(m, 0, y) for y in range(m))
        assert row == col == (0,) * (m - 1) + (1,)
        for shift in range(1, m):
            assert any(row[i] != row[(i + shift) % m] for i in range(m))
        for qx in range(3 * m + 1):
            for qy in range(3 * m + 1):
                if gcd(qx, qy) != 1:
                    continue
                codes, _ = return_cycle(m, qx % m, qy % m)
                residue = (qx + qy) % m
                expected = {0: (0,), 1: (0,) * (m - 2) + (1,),
                            m - 1: (0,) + (1,) * (m - 2)}.get(residue)
                if expected is None:
                    continue
                assert all(codes[j % len(codes)] == expected[j % len(expected)]
                           for j in range(4 * m))
                assert all(derivative(m, j * qx, j * qy) == expected[j % len(expected)]
                           for j in range(4 * m))
                n_rays += 1
        assert (derivative(m, 0, 0), derivative(m, 1, 0)) == (0, 0)
    print(f'INPUT: {n_blocks} block-periodicity checks; exact row/column period m for 4<=m<=12')
    print(f'RAYS: {n_rays} primitive directions in the three residue classes checked by literal enumeration')


def crt_pair(x, modulus, residue, next_modulus):
    assert gcd(modulus, next_modulus) == 1
    step = ((residue - x) * pow(modulus, -1, next_modulus)) % next_modulus
    combined = x + modulus * step
    assert combined % modulus == x % modulus
    assert combined % next_modulus == residue % next_modulus
    return combined, modulus * next_modulus


def obstruction(m, N):
    K = m * (m - 1) * factorial(N)
    cs, product_so_far = [], 1
    for _ in range(N):
        c = K * product_so_far - 1
        cs.append(c)
        product_so_far *= c
    Q = prod(cs)
    assert all(c > 1 and c % m == m - 1 and c % (m - 1) == m - 2 for c in cs)
    assert all(gcd(c, factorial(N)) == 1 for c in cs)
    assert all(gcd(c, d) == 1 for i, c in enumerate(cs) for d in cs[:i])
    a, modulus = (-Q) % m, m
    for j, c in enumerate(cs, 1):
        a, modulus = crt_pair(a, modulus, (-pow(j, -1, c)) % c, c)
    if a == 0:
        a += modulus
    assert a > 0 and Q > 0 and gcd(a, Q) == 1 and (a + Q) % m == 0
    pairs = []
    for ell, c in enumerate(cs, 1):
        assert (ell * a + 1) % c == 0
        x, y = ell * a + 1, ell * Q
        assert gcd(x, y) == gcd(x, Q) == c
        rx, ry = x // c, y // c
        assert gcd(rx, ry) == 1 and (rx + ry) % m == m - 1
        assert c % (m - 1) == m - 2
        pair = (derivative(m, ell * a, ell * Q), derivative(m, x, y))
        assert pair == (0, 1)
        pairs.append(pair)
    assert (0, 0) not in pairs
    return cs, a, Q


def family_checks():
    total = 0
    max_bits = 0
    composite_control = False
    for m in range(4, 13):
        for N in range(1, 8):
            cs, a, Q = obstruction(m, N)
            total += N
            max_bits = max(max_bits, a.bit_length(), Q.bit_length())
            if m == 4 and N == 2:
                assert cs[1] == 551 == 19 * 29
                composite_control = True
            print(f'FAMILY: m={m}, N={N}, exact-gcd checks={N}, adjacent derivative blocks=01^{N}')
    assert composite_control
    print(f'FAMILY TOTAL: 63 (m,N) pairs; {total} exact-gcd and adjacent-block checks; max coordinate bits={max_bits}')
    print('COMPOSITE CONTROL: m=4, N=2 uses c_2=551=19*29')


def main():
    source_fidelity()
    finite_periodicity_and_rays()
    family_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
