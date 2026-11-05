#!/usr/bin/env python3
"""Exact checks for Dietz--Flowers--Lockard Question 5.14.

A block (w, a) means a rows of length w, with distinct widths decreasing.
No numerical approximation or external dependency is used.
"""
from itertools import combinations


def transpose(blocks):
    """Ferrers reflection by constant-height column intervals."""
    total = 0
    columns = []
    for i, (width, count) in enumerate(blocks):
        assert width > 0 and count > 0
        lower = blocks[i + 1][0] if i + 1 < len(blocks) else 0
        assert lower < width
        total += count
        columns.append((total, width - lower))
    return tuple(reversed(columns))


def size(blocks):
    return sum(w * a for w, a in blocks)


def power_exponent(n, m):
    if n < 1:
        return None
    k = 0
    while n % m == 0:
        n //= m
        k += 1
    return k if n == 1 else None


def is_mary(blocks, m):
    return all(power_exponent(w, m) is not None for w, _ in blocks)


def blocks_from_rows(rows):
    return tuple((w, rows.count(w)) for w in sorted(set(rows), reverse=True))


def selfconj_from_widths(widths):
    """Candidate block diagram, checked by actual Ferrers transpose below."""
    assert tuple(sorted(set(widths))) == widths
    previous = 0
    counts = []
    for w in widths:
        counts.append(w - previous)
        previous = w
    return tuple(zip(reversed(widths), counts))


def check_congruence(m, blocks):
    assert blocks and transpose(blocks) == blocks and is_mary(blocks, m)
    L = blocks[0][0]
    K = power_exponent(L, m)
    a = [0] * (K + 1)
    for w, count in blocks:
        a[power_exponent(w, m)] = count
    assert sum(a) == L
    assert all(power_exponent(sum(a[j:]), m) is not None for j in range(K + 1))
    r = m - 1
    N = size(blocks)
    rhs = sum(a[j] * (m**j - 1) for j in range(K)) + (a[K] - 1) * (L - 1)
    assert N - (2 * L - 1) == rhs
    assert (N - (2 * L - 1)) % (r * r) == 0
    if m % 2:
        assert N % (2 * r) == 1
        e = power_exponent(N, m)
        assert e is None or e % 2 == 0
    return N


def small_partitions(n, bound):
    if n == 0:
        yield ()
    else:
        for first in range(min(n, bound), 0, -1):
            for rest in small_partitions(n - first, first):
                yield (first,) + rest


def main():
    # Definition 1.1's explicit example, including its conjugate and size.
    example = ((9, 1), (3, 2))
    assert size(example) == 15
    assert transpose(example) == ((3, 3), (1, 6))
    assert is_mary(example, 3) and is_mary(transpose(example), 3)
    print('FIDELITY Definition 1.1: (1,2,0)_3 -> (3,6)_3; size=15')

    # Proposition 1.3: compare direct transpose with multiplicity partial sums
    # for every nonzero tuple in {0,...,4}^5, for each base 2,...,9.
    from itertools import product
    fidelity_cases = 0
    for m in range(2, 10):
        for a in product(range(5), repeat=5):
            if not any(a):
                continue
            K = max(j for j, count in enumerate(a) if count)
            blocks = tuple((m**j, a[j]) for j in range(K, -1, -1) if a[j])
            lhs = is_mary(transpose(blocks), m)
            rhs = all(power_exponent(sum(a[j:K + 1]), m) is not None for j in range(K + 1))
            assert lhs == rhs
            assert transpose(transpose(blocks)) == blocks
            fidelity_cases += 1
    print(f'FIDELITY Proposition 1.3: {fidelity_cases} exact multiplicity tuples OK')

    # Independent direct cell reflection for every partition of n <= 35.
    partitions = selfconj = mary_selfconj = 0
    for n in range(1, 36):
        for rows in small_partitions(n, n):
            partitions += 1
            blocks = blocks_from_rows(rows)
            cells = {(i, j) for i, w in enumerate(rows) for j in range(w)}
            reflected = {(j, i) for i, j in cells}
            tb = transpose(blocks)
            tr = tuple(w for w, a in tb for _ in range(a))
            assert reflected == {(i, j) for i, w in enumerate(tr) for j in range(w)}
            if reflected != cells:
                continue
            selfconj += 1
            for m in range(2, 36):
                if is_mary(blocks, m):
                    check_congruence(m, blocks)
                    mary_selfconj += 1
    print(f'DIRECT n<=35: partitions={partitions}; selfconjugate={selfconj}; m-ary pairs={mary_selfconj}; obstruction OK')

    # Enumerate all 511 nonempty sets of widths from {m^0,...,m^8}.
    # Each constructed diagram is checked by transpose, not assumed symmetric.
    odd_count = all_count = 0
    for m in range(2, 23):
        powers = tuple(m**j for j in range(9))
        count = 0
        for k in range(1, 10):
            for widths in combinations(powers, k):
                blocks = selfconj_from_widths(widths)
                check_congruence(m, blocks)
                count += 1
        assert count == 511
        all_count += count
        if m in range(3, 22, 2):
            odd_count += count
            print(f'ODD base={m}: selfconjugate diagrams={count}; square-modulus and parity obstruction OK')
    assert odd_count == 5110 and all_count == 10731
    print(f'COMPRESSED: total={all_count}; designated odd-base cases={odd_count}')

    # Published square family and the necessity of restricting to odd bases.
    square_cases = 0
    for m in range(2, 23):
        for t in range(9):
            blocks = ((m**t, m**t),)
            assert transpose(blocks) == blocks and size(blocks) == m**(2 * t)
            check_congruence(m, blocks)
            square_cases += 1
    even_example = ((4, 1), (2, 1), (1, 2))
    assert transpose(even_example) == even_example
    assert is_mary(even_example, 2) and size(even_example) == 2**3
    print(f'FIDELITY Discussion 5.11 square family: {square_cases} cases; even-base boundary (4,2,1,1) has size 8')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
