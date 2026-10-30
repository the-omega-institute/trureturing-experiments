#!/usr/bin/env python3
"""Exact finite checks for Kurkov's A232559 conjecture; not a proof.

Published fidelity data: OEIS A232559 r63, A059894 r97, A348366 r23,
A072649 r104, A130312 r39. All arithmetic is integral. No network required.
"""
from collections import deque

LIMIT = 50000
MAX_DEPTH = 21


def fifo(count):
    """Literal source definition, including odd increments and immediate deduplication."""
    queue, seen, values, depths = deque([(1, 0)]), {1}, [], []
    while len(values) < count:
        x, depth = queue.popleft()
        values.append(x)
        depths.append(depth)
        for y in (x + 1, 2 * x):
            if y not in seen:
                seen.add(y)
                queue.append((y, depth + 1))
    return values, depths


def reverse_complement(n):
    assert n > 0
    suffix = bin(n)[3:]
    return int('1' + ''.join('1' if b == '0' else '0' for b in reversed(suffix)), 2)


def fibonacci_data(count):
    fib = [0, 1]
    while fib[-1] <= count:
        fib.append(fib[-1] + fib[-2])
    h, q, level = [0], [0], 1
    for n in range(1, count + 1):
        while fib[level + 2] <= n:
            level += 1
        assert fib[level + 1] <= n < fib[level + 2]
        h.append(level)
        q.append(fib[level])
    return fib, h, q


def recurrence(count, h, q):
    b = [0, 1]
    for n in range(2, count + 1):
        m = n - q[n]
        assert 1 <= m < n
        assert h[m] < h[n]
        b.append(2 * b[m] + (h[n] - h[m]) % 2)
    return b


def weighted_words(max_depth):
    """Enumerate all suffixes of cost <= max_depth; sort actual operation paths.

    This does not use the proposed row concatenation recurrence. D has cost 1,
    and a suffix 1 is encoded by DA, with A preceding D in sibling order.
    """
    rows = [[] for _ in range(max_depth + 1)]
    stack = [('', 0, ())]
    while stack:
        word, cost, path = stack.pop()
        rows[cost].append((path, word))
        if cost + 1 <= max_depth:
            stack.append((word + '0', cost + 1, path + (1,)))
        if cost + 2 <= max_depth:
            stack.append((word + '1', cost + 2, path + (1, 0)))
    return [[word for _, word in sorted(row)] for row in rows]


def main():
    a, depths = fifo(LIMIT)
    fib, h, q = fibonacci_data(LIMIT)
    b = recurrence(LIMIT, h, q)
    computed = {'232559': a, '059894': [reverse_complement(n) for n in range(1, len(PUBLISHED['059894']) + 1)],
                '348366': b, '072649': h[1:], '130312': q[1:]}
    for name, published in PUBLISHED.items():
        assert computed[name][:len(published)] == published, name
        print(f'FIDELITY A{name}: {len(published)} published terms matched')
    source_rows = [[1], [2], [3, 4], [6, 5, 8], [7, 12, 10, 9, 16]]
    for d, row in enumerate(source_rows):
        assert [x for x, depth in zip(a[:12], depths[:12]) if depth == d] == row
    print('FIDELITY A232559 COMMENTS: first 5 ordered generations matched')
    for n, x in enumerate(a, 1):
        assert x == reverse_complement(b[n]), (n, x, b[n])
        assert reverse_complement(reverse_complement(x)) == x
        assert depths[n - 1] == x.bit_length() + x.bit_count() - 2
        assert h[n] == depths[n - 1] + 1
    print(f'IDENTITY: all indices 1..{LIMIT} matched; involution and depth also checked')
    rows = weighted_words(MAX_DEPTH)
    flattened = []
    for d, words in enumerate(rows):
        assert len(words) == fib[d + 1]
        values = [int('1' + w, 2) for w in words]
        assert values == a[fib[d + 2] - 1:fib[d + 3] - 1]
        assert all(len(w) + w.count('1') == d for w in words)
        if d >= 2:
            assert words == ['1' + w for w in rows[d - 2]] + ['0' + w for w in rows[d - 1]]
            conjugated = [reverse_complement(x) for x in values]
            expected = [2 * reverse_complement(int('1' + w, 2)) for w in rows[d - 2]]
            expected += [2 * reverse_complement(int('1' + w, 2)) + 1 for w in rows[d - 1]]
            assert conjugated == expected
        flattened.extend(values)
    assert len(flattened) == fib[MAX_DEPTH + 3] - 1
    assert flattened == a[:len(flattened)]
    print(f'ORDERED ROWS: depths 0..{MAX_DEPTH}, {len(flattened)} nodes; exact order, counts, offsets, and conjugation matched')
    print('ALL_CHECKS_OK')


PUBLISHED = {'232559': [1, 2, 3, 4, 6, 5, 8, 7, 12, 10, 9, 16, 14, 13, 24, 11, 20, 18, 17, 32, 15, 28, 26, 25, 48, 22, 21, 40, 19, 36, 34, 33, 64, 30, 29, 56, 27, 52, 50, 49, 96, 23, 44, 42, 41, 80, 38, 37, 72, 35, 68, 66, 65, 128, 31, 60, 58, 57, 112, 54, 53, 104, 51, 100, 98, 97], '059894': [1, 3, 2, 7, 5, 6, 4, 15, 11, 13, 9, 14, 10, 12, 8, 31, 23, 27, 19, 29, 21, 25, 17, 30, 22, 26, 18, 28, 20, 24, 16, 63, 47, 55, 39, 59, 43, 51, 35, 61, 45, 53, 37, 57, 41, 49, 33, 62, 46, 54, 38, 58, 42, 50, 34, 60, 44, 52, 36, 56, 40, 48, 32, 127, 95, 111, 79, 119, 87, 103, 71], '348366': [0, 1, 3, 2, 7, 6, 5, 15, 4, 14, 13, 11, 31, 12, 10, 30, 9, 29, 27, 23, 63, 8, 28, 26, 22, 62, 25, 21, 61, 19, 59, 55, 47, 127, 24, 20, 60, 18, 58, 54, 46, 126, 17, 57, 53, 45, 125, 51, 43, 123, 39, 119, 111, 95, 255, 16, 56, 52, 44, 124, 50, 42, 122, 38, 118], '072649': [1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10], '130312': [1, 1, 2, 2, 3, 3, 3, 5, 5, 5, 5, 5, 8, 8, 8, 8, 8, 8, 8, 8, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34]}

if __name__ == '__main__':
    main()

