#!/usr/bin/env python3
"""Exact checks for Barker's two OEIS A268943 conjectures; no Lean proof.

The finite objects are words, not sequences defined by the conjectured recurrence.
Run beside oeis-b268943.txt using Python 3; the standard library suffices.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import product
from pathlib import Path


def event_labels(word):
    return [right for left, right in zip(word, word[1:]) if left == right]


def legal(word, q, *, unequal=False):
    events = event_labels(word)
    if unequal:
        return all(b != a for a, b in zip(events, events[1:]))
    return all(b == (a + 1) % q for a, b in zip(events, events[1:]))


def literal_state_dp(q, limit):
    """Retain both actual last letter and last equality-event label (or None)."""
    states = {(a, None): 1 for a in range(q)}
    totals, partitions = [1], [None]
    for n in range(1, limit + 1):
        totals.append(sum(states.values()))
        u = sum(v for (last, event), v in states.items() if event is None)
        a = sum(v for (last, event), v in states.items()
                if event is not None and last == (event + 1) % q)
        b = sum(v for (last, event), v in states.items()
                if event is not None and last != (event + 1) % q)
        partitions.append((u, a, b))
        nxt = defaultdict(int)
        for (last, event), count in states.items():
            for letter in range(q):
                if letter == last:
                    if event is None or letter == (event + 1) % q:
                        nxt[(letter, letter)] += count
                else:
                    nxt[(letter, event)] += count
        states = nxt
    return totals, partitions


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def ratio_coefficients(num, den, limit):
    assert den[0] == 1
    out = []
    for n in range(limit + 1):
        value = num[n] if n < len(num) else 0
        value -= sum(den[i] * out[n - i]
                     for i in range(1, min(n, len(den) - 1) + 1))
        out.append(value)
    return out


def qmul(x, y):
    """Multiplication in Q[s]/(s^2 - 17), with unique pair representation."""
    a, b = x
    c, d = y
    return a * c + 17 * b * d, a * d + b * c


def radical_values(limit):
    plus = minus = (Fraction(1), Fraction(0))
    out = []
    for n in range(limit + 1):
        x = qmul((Fraction(153), Fraction(-37)), minus)
        y = qmul((Fraction(153), Fraction(37)), plus)
        out.append(((x[0] + y[0] - 272 * 7**n) / 34,
                    (x[1] + y[1]) / 34))
        plus = qmul(plus, (3, 1))
        minus = qmul(minus, (3, -1))
    return out


def main():
    limit = 210
    all_totals = {}
    comparisons = 0
    transition_checks = 0
    for q in range(2, 11):
        counts, parts = literal_state_dp(q, limit)
        den = poly_mul([1, -(q - 1)], [1, -(q - 2), -q])
        num = [0, q, -q * (q - 3), -q * (q - 1)]
        gf = ratio_coefficients(num, den, limit)
        assert counts[0] == 1 and gf[0] == 0
        assert counts[1:] == gf[1:], q
        comparisons += len(counts)
        assert parts[1] == (q, 0, 0)
        for n in range(1, limit):
            u, a, b = parts[n]
            assert parts[n + 1] == ((q - 1) * u, b,
                                     u + q * a + (q - 2) * b), (q, n)
            transition_checks += 1
        all_totals[q] = counts
    print(f'Literal full-state DP vs rational GF: q=2..10, n=0..{limit}, '
          f'{comparisons} comparisons (at n=0 use 1+GF): PASS')
    print(f'Three-state transition vs full-state DP: {transition_checks} checks: PASS')

    exhaustive_words = exhaustive_counts = 0
    bounds = {2: 10, 3: 8, 4: 7, 5: 7, 8: 6}
    for q, bound in bounds.items():
        for n in range(bound + 1):
            actual = sum(legal(word, q) for word in product(range(q), repeat=n))
            assert actual == all_totals[q][n], (q, n, actual)
            exhaustive_words += q**n
            exhaustive_counts += 1
    print(f'Exhaustive words: {exhaustive_words}; {exhaustive_counts} counts; '
          f'bounds {bounds}: PASS')

    expected_prefix = [8, 64, 504, 3928, 30344, 232696, 1773384, 13443064,
                       101433800, 762265720, 5707893576, 42605289208,
                       317113497800, 2354253598072, 17437541654088,
                       128885063291896, 950791205063624, 7001691181273720,
                       51477520840048968]
    counts = all_totals[8]
    assert counts[1:20] == expected_prefix
    bfile = Path(__file__).with_name('oeis-b268943.txt')
    source_terms = {}
    for line in bfile.read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            n, value = map(int, line.split())
            source_terms[n] = value
    assert sorted(source_terms) == list(range(1, 211))
    assert all(counts[n] == value for n, value in source_terms.items())
    print('Source fidelity: 19 displayed terms and all 210 b-file terms: PASS')

    rows = [
        [0, 7, 3, 6, 2, 5, 7, 0, 1, 2, 2, 0, 6, 3, 5, 7],
        [7, 6, 5, 7, 1, 2, 7, 5, 1, 1, 4, 5, 4, 3, 5, 5],
        [2, 3, 5, 6, 0, 0, 3, 5, 2, 3, 7, 2, 3, 1, 2, 6],
        [1, 2, 0, 5, 3, 7, 5, 1, 0, 5, 5, 4, 5, 5, 5, 1],
        [1, 0, 3, 3, 3, 0, 6, 5, 7, 1, 3, 4, 4, 6, 4, 6],
    ]
    assert all(legal(word, 8) for word in zip(*rows))
    assert legal((6, 4, 6, 5, 4), 8)
    assert legal((7, 7, 0, 0), 8)  # modulo-eight wrap
    assert not legal((0, 0, 0), 8)  # overlapping equality events
    assert not legal((0, 0, 2, 2), 8)
    print('Source fidelity: all 16 example columns; wrap and overlapping-event guards: PASS')

    num, den = [0, 8, -40, -56], [1, -13, 34, 56]
    assert ratio_coefficients(num, den, limit)[1:] == counts[1:]
    assert ratio_coefficients([1, -5, -6], den, limit) == counts
    assert poly_mul([-7, 1], [-8, -6, 1]) == [56, 34, -13, 1]
    for n, value in enumerate(radical_values(limit)):
        assert value == (Fraction(counts[n]), Fraction(0)), (n, value)
    print('Barker radical expression in Q[s]/(s^2-17), n=0..210: PASS')
    print('Characteristic polynomial factorization and 1+GF coefficients: PASS')
    unequal = sum(legal(w, 8, unequal=True) for w in product(range(8), repeat=4))
    assert unequal == 3976 and counts[4] == 3928
    print(f'Fried unequal-event relation vs cyclic-successor relation at q=8,n=4: '
          f'{unequal} vs {counts[4]}: DISTINCT')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
