#!/usr/bin/env python3
"""Exact finite checks for arXiv:2610.00765v1, Section 7, Question (1).

Standard library only. No numerical approximation and no floating point.
Finite tests do not prove the unbounded proposed settlement.
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import gcd, isqrt, prod


def factors(n):
    answer = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            answer[p] = answer.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        answer[n] = 1
    return answer


def integer_square(q):
    return q.denominator == 1 and q.numerator >= 0 and isqrt(q.numerator) ** 2 == q.numerator


def actual_sequence(h, start, last):
    x = [Fraction(start)]
    total = x[0]
    for _ in range(last):
        x.append(Fraction(h * h, 1) / total)
        total += x[-1]
    return x


def actual_failures(h, start, last=7):
    x = actual_sequence(h, start, last)
    values = {n: Fraction(2, h) * x[n].numerator + x[n].denominator
              for n in range(2, last + 1)}
    return {n for n, value in values.items() if not integer_square(value)}, values


def hit_times(p):
    """Independent forward field walks; 0 denotes infinity/no positive hit."""
    answer = {}
    for z in range(1, p):
        seen = set()
        w, t = z, 0
        while w and w not in seen:
            seen.add(w)
            w = (w + pow(w, -1, p)) % p
            t += 1
        answer[z] = t if w == 0 else 0
        assert answer[z] < p
    return answer


def backward_sets(p):
    """The proposed quadratic-preimage construction, independent of walks."""
    previous, out = {0}, {}
    for t in range(1, p + 1):
        current = {z for z in range(1, p)
                   if any((z * z + 1 - w * z) % p == 0 for w in previous)}
        if current:
            out[t] = current
        previous = current
    assert not previous
    return out


def local_data(h):
    out = []
    for p, e in factors(h).items():
        if p == 2 or p % 4 == 1:
            times = hit_times(p)
            backward = backward_sets(p)
            assert {t: {z for z in times if times[z] == t}
                    for t in sorted(set(times.values()) - {0})} == backward
            m = p ** (e + 1)
            weights = Counter(times[z] for z in times if times[z])
            weights[0] = m - sum(weights.values())
            residue_labels = {p ** e * (h // p ** e) * z % m: t
                              for z, t in times.items() if t}
            assert len(residue_labels) == sum(v for t, v in weights.items() if t)
            out.append((p, m, weights, residue_labels))
    return out


def proposed(h, start, data):
    return {labels.get(start % m, 0) + 1 for _, m, _, labels in data
            if labels.get(start % m, 0)}


def weighted_formula(data):
    counts = Counter()
    for choices in product(*(list(weights.items()) for _, _, weights, _ in data)):
        labels = {t for t, _ in choices} - {0}
        counts[len(labels)] += prod(w for _, w in choices)
    return counts


def homogeneous_orbit(h, start):
    """Full finite-state orbit modulo h, from the *actual* reduced start.

    Detect the eventual repeated state, and reject a cycle containing any
    obstruction. This checks all positive times in this modular orbit; it
    does not assume the proposed prime-power labels or their CRT product.
    """
    g = gcd(h, start)
    u, v = start // g % h, h // g % h
    seen, events, n = {}, [], 0
    while (u, v) not in seen:
        seen[u, v] = n
        if n >= 1 and gcd(h, u) > 1:
            events.append(n + 1)
        u, v = (u * u + v * v) % h, u * v % h
        n += 1
    cycle_start = seen[u, v]
    assert all(event - 1 < cycle_start for event in events), (h, start, cycle_start, events)
    return set(events)


def main():
    # Source fidelity: Example 2.4, and Theorem B's h=x0=10 example.
    x = actual_sequence(3, 1, 4)
    assert x == list(map(Fraction, [1, 9])) + [Fraction(9, 10), Fraction(90, 109), Fraction(9810, 12781)]
    bad, values = actual_failures(3, 1, 4)
    assert not bad and list(values.values()) == [16, 169, 19321]
    bad10, values10 = actual_failures(10, 10, 5)
    assert bad10 == {2, 3}
    assert list(values10.values()) == [Fraction(2), Fraction(9, 5), Fraction(49), Fraction(1521)]
    print('FIDELITY h=3 x0=1: x=[1,9,9/10,90/109,9810/12781]; square readouts=[16,169,19321]')
    print('FIDELITY h=10 x0=10: readouts=[2,9/5,49,1521]; failures=[2,3]')

    # Direct unreduced recurrence vs Fraction's intrinsic reduced coordinates.
    direct_cases = 0
    for h in range(1, 81):
        data = local_data(h)
        for start in range(1, 161):
            bad, values = actual_failures(h, start)
            predicted = proposed(h, start, data)
            assert bad == {n for n in predicted if n <= 7}, (h, start, bad, predicted)
            g = gcd(h, start)
            u, v = start // g, h // g
            for t in range(1, 7):
                old_u, old_v = u, v
                u, v = u * u + v * v, u * v
                assert gcd(u, v) == 1
                assert values[t + 1] == Fraction((old_u + old_v) ** 2, gcd(h, u))
                assert integer_square(values[t + 1]) == (gcd(h, u) == 1)
            direct_cases += 1
    print(f'DIRECT exact Fraction checks: h=1..80, x0=1..160, n=2..7; {direct_cases} starts')

    class_cases = 0
    for h in list(range(1, 41)) + [50, 65]:
        data = local_data(h)
        M = prod(m for _, m, _, _ in data)
        formula = weighted_formula(data)
        assert sum(formula.values()) == M
        observed, representatives = Counter(), {}
        for start in range(1, M + 1):
            predicted = proposed(h, start, data)
            assert homogeneous_orbit(h, start) == predicted, (h, start)
            assert homogeneous_orbit(h, start + M) == predicted, (h, start + M)
            observed[len(predicted)] += 1
            representatives.setdefault(len(predicted), []).append(start)
        assert observed == formula, (h, observed, formula)
        for N in [0, 1, M - 1, M, M + 1, 2 * M + M // 2, 3 * M + 7]:
            for j in range(len(data) + 2):
                brute = sum(len(proposed(h, start, data)) == j for start in range(1, N + 1))
                q, s = divmod(N, M)
                exact = q * formula[j] + sum(r <= s for r in representatives.get(j, []))
                assert brute == exact
                assert abs(Fraction(brute) - Fraction(N * formula[j], M)) <= M
        class_cases += M
    print(f'FINITE-STATE/CRT checks: every class for h=1..40,50,65, and its shift by M; {class_cases} classes')

    h, data = 65, local_data(65)
    counts = weighted_formula(data)
    two = [r for r in range(1, 4226) if len(proposed(h, r, data)) == 2]
    collisions = [r for r in range(1, 4226)
                  if sum(r % m in labels for _, m, _, labels in data) == 2
                  and len(proposed(h, r, data)) == 1]
    assert (counts[0], counts[1], counts[2]) == (3507, 714, 4)
    assert two == [1365, 2015, 2210, 2860]
    assert collisions == [520, 1170, 3055, 3705]
    bad, values = actual_failures(65, 1365, 7)
    assert bad == {2, 3}
    assert values[2] == Fraction(484, 13) and values[3] == Fraction(214369, 5)
    print('h=65: M=4225; counts=(3507,714,4); densities=(3507/4225,714/4225,4/4225)')
    print('h=65 two-failure classes:', two)
    print('h=65 coincident-prime classes:', collisions)
    print('h=65 x0=1365: failing readouts at n=2,3 are 484/13 and 214369/5')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
