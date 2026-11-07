#!/usr/bin/env python3
"""Exact, independent width-layer and divisor-block checks for OEIS A401249.

The geometric side uses A237048's triangular-offset columns and A249223's
alternating partial sums, reflected as prescribed by A262045. It does not
use a divisor block, or an odd-divisor count, to define a geometric part.
Only Python's integer arithmetic is used. No third-party modules are needed.
"""

from itertools import groupby
from math import isqrt


def divisors(n):
    low = [d for d in range(1, isqrt(n) + 1) if n % d == 0]
    return low + [n // d for d in reversed(low) if d * d != n]


def blocks(n):
    result = []
    for d in divisors(n):
        if not result or d > 2 * result[-1][-1]:
            result.append([])
        result[-1].append(d)
    return result


def half_widths(n):
    """Literal A237048 definition, then A249223 alternating partial sums."""
    r = (isqrt(8 * n + 1) - 1) // 2
    result = []
    width = 0
    for k in range(1, r + 1):
        offset = k * (k + 1) // 2
        b = int(n >= offset and (n - offset) % k == 0)
        width += b if k % 2 else -b
        assert width >= 0, (n, k, width)
        result.append(width)
    return result


def widths(n):
    h = half_widths(n)
    return h + h[::-1]


def positive_parts(w):
    return [list(v) for positive, v in groupby(w, key=lambda x: x > 0)
            if positive]


def components_at_level(w, level):
    return sum(1 for positive, _ in groupby(w, key=lambda x: x >= level)
               if positive)


def layer_counts(w):
    return [components_at_level(w, k) for k in range(1, max(w) + 1)]


def part_counts(n):
    return [sum(layer_counts(p)) for p in positive_parts(widths(n))]


def odd_counts(n):
    return [sum(d % 2 for d in b) for b in blocks(n)]


def main():
    # Fidelity: the primary source's fully worked example.
    assert divisors(70) == [1, 2, 5, 7, 10, 14, 35, 70]
    assert blocks(70) == [[1, 2], [5, 7, 10, 14], [35, 70]]
    assert odd_counts(70) == [1, 2, 1]
    assert min(odd_counts(70)) == 1
    print('FIDELITY A401249 n=70: blocks=[[1, 2], [5, 7, 10, 14], [35, 70]]; odd counts=[1, 2, 1]; minimum=1')

    # Fidelity of the geometry: examples in A262045 and A279387.
    assert widths(45) == [1, 0, 1, 1, 2, 1, 1, 1, 2,
                          2, 1, 1, 1, 2, 1, 1, 0, 1]
    examples = {12: [1, 1], 15: [3, 1], 18: [1, 2],
                45: [3, 3], 90: [1, 3, 2], 360: [1, 1, 1, 1, 2]}
    for n, expected in examples.items():
        assert layer_counts(widths(n)) == expected, (n, expected)
        print(f'FIDELITY A279387 n={n}: subparts by layer={expected}; subparts by part={part_counts(n)}')
    assert part_counts(12) == [2]
    assert part_counts(15) == [1, 2, 1]  # A280940's worked examples.

    # Primary source revision #22's entire published DATA (105 terms).
    data = [
        1,1,1,1,1,2,1,1,1,1,1,2,1,1,1,1,1,3,1,2,1,1,1,2,1,1,1,2,1,4,1,1,1,1,
        1,3,1,1,1,2,1,4,1,1,1,1,1,2,1,1,1,1,1,4,1,2,1,1,1,4,1,1,1,1,1,4,1,1,
        1,1,1,3,1,1,1,1,1,2,1,2,1,1,1,4,1,1,1,2,1,6,1,1,1,1,1,2,1,1,1,3,1,2,1,2,1,
    ]
    for n, expected in enumerate(data, 1):
        assert min(odd_counts(n)) == expected, (n, expected)
        assert min(part_counts(n)) == expected, (n, expected)
    print(f'FIDELITY A401249 revision #22: all {len(data)} DATA terms reproduced by both definitions')

    limit = 10000
    totals = dict(parts=0, subparts=0, half_width_entries=0,
                  central=0, noncentral=0)
    for n in range(1, limit + 1):
        bs = blocks(n)
        h = half_widths(n)
        w = h + h[::-1]
        ps = positive_parts(w)
        lhs = [sum(layer_counts(p)) for p in ps]
        rhs = [sum(d % 2 for d in b) for b in bs]
        assert lhs and rhs
        assert lhs == rhs, ('pointwise', n, lhs, rhs)
        assert min(lhs) == min(rhs), ('minimum', n, lhs, rhs)
        assert sum(lhs) == sum(d % 2 for d in divisors(n)), ('total', n)
        # An independently stated known identity, Hoft A249223 (2024).
        ds = divisors(n)
        assert h == [sum(k < 2 * d and d <= k for d in ds)
                     for k in range(1, len(h) + 1)], ('window', n)
        # Layer components equal positive variation, including each part's
        # initial rise from zero. This checks the proposed local-count lemma.
        for p, count in zip(ps, lhs):
            assert count == sum(max(0, b - a) for a, b in zip([0] + p, p))
        # Every complete dyadic chain lies in one actual block.
        m, q = 0, n
        while q % 2 == 0:
            m, q = m + 1, q // 2
        for b in bs:
            assert len(b) == (m + 1) * sum(d % 2 for d in b)
            for u in b:
                if u % 2:
                    assert all(2**j * u in b for j in range(m + 1))
        totals['parts'] += len(lhs)
        totals['subparts'] += sum(lhs)
        totals['half_width_entries'] += len(h)
        totals['central'] += int(h[-1] > 0)
        totals['noncentral'] += len(lhs) - int(h[-1] > 0)
    print(f'EXACT CHECK n=1..{limit}: ordered per-part counts = ordered per-block odd counts; minima equal')
    print('EXACT CHECK: total subparts = odd divisors; alternating widths = half-open divisor-window widths')
    print('EXACT CHECK: layer-component counts = positive variation; complete dyadic chains stay in one block')
    print('COVERAGE ' + ' '.join(f'{k}={v}' for k, v in totals.items()))
    print('SCOPE: finite exact evidence only; no universal proof or Lean probe was performed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
