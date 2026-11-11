#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
"""Exact finite checks of Bradford's formula and the four-generator construction.

These checks use finite truncations of the sparse lamps. They do not prove the
infinite-family statement. All matrix entries and coordinates are integers.
"""

from itertools import product


ONE = (1, 0, 0, 1)


def matrices(modulus):
    def mul(x, y):
        a, b, c, d = x
        e, f, g, h = y
        return ((a * e + b * g) % modulus, (a * f + b * h) % modulus,
                (c * e + d * g) % modulus, (c * f + d * h) % modulus)

    def inv(x):
        a, b, c, d = x
        assert (a * d - b * c) % modulus == 1
        return (d, -b % modulus, -c % modulus, a)

    return mul, inv


def words(rank, radius):
    alphabet = tuple(range(1, rank + 1)) + tuple(range(-rank, 0))

    def extend(prefix):
        if prefix:
            yield prefix
        if len(prefix) < radius:
            for letter in alphabet:
                if not prefix or letter != -prefix[-1]:
                    yield from extend(prefix + (letter,))

    yield from extend(())


def reduce_word(word):
    stack = []
    for letter in word:
        if stack and stack[-1] == -letter:
            stack.pop()
        else:
            stack.append(letter)
    return tuple(stack)


def inverse_word(word):
    return tuple(-a for a in reversed(word))


def z_word(j):
    return (1,) * j + (2,) + (-1,) * j + (-2,)


def evaluate(word, generators, mul, inv, identity):
    value = identity
    for letter in word:
        factor = generators[abs(letter) - 1]
        value = mul(value, factor if letter > 0 else inv(factor))
    return value


def positions(size):
    p, q = [0], [0]
    while len(p) < size:
        p.append(p[-1] + q[-1] + 1)
        q.append(p[-1] + q[-1] + 1)
    return p, q


def check_sparse_positions():
    p, q = positions(16)
    mixed = [q[j] - p[k] for j, k in product(range(16), repeat=2)]
    assert len(mixed) == len(set(mixed))
    for s in (p, q):
        differences = [s[j] - s[k] for j, k in product(range(16), repeat=2) if j != k]
        assert len(differences) == len(set(differences))
    print("FIDELITY sparse-position lemma: 16 indices, 256 distinct mixed differences")


def lamp_operations(modulus):
    mmul, minv = matrices(modulus)

    def clean(f):
        return {x: v for x, v in f.items() if v != ONE}

    def shift(f, k):
        return {x + k: v for x, v in f.items()}

    def pmul(f, g):
        return clean({x: mmul(f.get(x, ONE), g.get(x, ONE)) for x in f.keys() | g.keys()})

    def pinv(f):
        return clean({x: minv(v) for x, v in f.items()})

    def wmul(a, b):
        f, k = a
        g, l = b
        return pmul(f, shift(g, k)), k + l

    def winv(a):
        f, k = a
        return shift(pinv(f), -k), -k

    return mmul, minv, shift, pmul, pinv, wmul, winv


def check_source_formula():
    # Sanov's integral matrices, reduced modulo a finite integer modulus.
    # Only the explicitly enumerated words are claimed as checked here.
    ell = 6
    modulus = 2 * 3**ell + 1
    a, b = (1, 2, 0, 1), (1, 0, 2, 1)
    mmul, minv, shift, _, _, wmul, winv = lamp_operations(modulus)
    p, q = positions(8)
    g = {x: a for x in p}
    h = {x: b for x in q}
    index = 4
    aligned = (shift(g, q[index] - p[index]), 0)
    count = 0
    for word in words(2, ell):
        direct = evaluate(word, (a, b), mmul, minv, ONE)
        assert direct != ONE
        lifted = evaluate(word, (aligned, (h, 0)), wmul, winv, ({}, 0))
        assert lifted[0].get(q[index], ONE) == direct
        count += 1
    assert count == 2 * (3**ell - 1)
    print(f"FIDELITY Bradford evaluation formula: {count} reduced words of lengths 1..{ell}, modulus {modulus}")


def check_radius(radius):
    M = 2 * radius + 1
    m = 3 * M
    c = -radius - 1
    ell = radius * (2 * m + 2)
    modulus = 2 * 3**ell + 1
    a, b = (1, 2, 0, 1), (1, 0, 2, 1)
    mmul, minv, shift, pmul, pinv, wmul, winv = lamp_operations(modulus)
    p, q = positions(12)
    g, h = {x: a for x in p}, {x: b for x in q}
    index = 5
    aligned_g, aligned_h = shift(g, -p[index]), shift(h, -q[index])
    z = {}
    for j in range(1, m + 1):
        word = z_word(j)
        assert len(word) <= 2 * m + 2
        assert sum(x == 1 for x in word) == sum(x == -1 for x in word)
        assert sum(x == 2 for x in word) == sum(x == -2 for x in word)
        value = evaluate(word, (a, b), mmul, minv, ONE)
        assert value != ONE
        single = evaluate(word, ((aligned_g, 0), (aligned_h, 0)), wmul, winv, ({}, 0))
        assert single == ({0: value}, 0)
        z[j] = value

    prescribed = {}
    constructed = {}
    for i, k in product(range(3), range(-radius, radius + 1)):
        site = c - k - i * M
        j = 1 + i * M + k + radius
        assert site < 0 and site not in prescribed
        prescribed[site] = z[j]
        # The factor is a translate of a word in aligned original generators.
        single = evaluate(z_word(j), ((aligned_g, 0), (aligned_h, 0)), wmul, winv, ({}, 0))
        constructed = pmul(constructed, shift(single[0], site))
    assert constructed == prescribed
    y = constructed
    generators = (({}, 1), (y, 0), (pmul(g, shift(y, M)), 0), (pmul(h, shift(y, 2 * M)), 0))
    assert pmul(generators[2][0], pinv(shift(y, M))) == g
    assert pmul(generators[3][0], pinv(shift(y, 2 * M))) == h

    count = zero_shift = 0
    for word in words(4, radius):
        actual = evaluate(word, generators, wmul, winv, ({}, 0))
        shift_sum, labels, expansion = 0, [], ()
        for letter in word:
            if abs(letter) == 1:
                shift_sum += 1 if letter > 0 else -1
            else:
                i = abs(letter) - 2
                j = 1 + i * M + shift_sum + radius
                assert -radius <= shift_sum <= radius and 1 <= j <= m
                labels.append(j if letter > 0 else -j)
                expansion += z_word(j) if letter > 0 else inverse_word(z_word(j))
        assert actual[1] == shift_sum
        assert actual != ({}, 0)
        if shift_sum == 0:
            zero_shift += 1
            assert labels and all(x != -y for x, y in zip(labels, labels[1:]))
            assert reduce_word(expansion)
            assert len(expansion) <= ell
            expected = evaluate(expansion, (a, b), mmul, minv, ONE)
            assert actual[0].get(c, ONE) == expected != ONE
        count += 1
    assert count == 8 * (7**radius - 1) // 6
    print(f"CONSTRUCTION r={radius}: M={M}, m={m}, ell={ell}, all {count} reduced words pass; {zero_shift} have zero shift")
    print(f"MEMBERSHIP r={radius}: {m} single-site words and translated factors recover original generators")


def main():
    check_sparse_positions()
    check_source_formula()
    for radius in (2, 3):
        check_radius(radius)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
