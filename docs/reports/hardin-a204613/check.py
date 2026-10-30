#!/usr/bin/env python3
"""Exact finite check for OEIS A204613 (Hardin, revision 7)."""
from collections import Counter
from itertools import product
from pathlib import Path

COEFFICIENTS = [
    36, -406, 404, 23747, -145816, -129968, 3736568,
    -7902612, -29701104, 128249776, 34934976, -780114560,
    607369728, 2210927104, -3052507136, -2863537152,
    5709414400, 1514283008, -4781309952, -281149440,
    1658847232, 94371840, -150994944,
]

ROWS = list(product((0, 1), repeat=6))


def determinant_vector(a, b):
    return tuple(a[j] * b[j + 1] - a[j + 1] * b[j] for j in range(5))


def monotone(d):
    return all(d[j] <= d[j + 1] for j in range(4))


def count_sequence(length):
    """Return a(1), ..., a(length) by exact weighted state DP."""
    weights = Counter(
        (bottom, determinant_vector(top, bottom))
        for top in ROWS
        for bottom in ROWS
        if monotone(determinant_vector(top, bottom))
    )
    states = sorted(weights)
    state_index = {state: i for i, state in enumerate(states)}
    edges = []
    for bottom, d in states:
        targets = []
        for next_bottom in ROWS:
            next_d = determinant_vector(bottom, next_bottom)
            if monotone(next_d) and all(x <= y for x, y in zip(d, next_d)):
                targets.append(state_index[(next_bottom, next_d)])
        edges.append(targets)

    values = [weights[state] for state in states]
    sequence = [sum(values)]
    for _ in range(length - 1):
        next_values = [0] * len(states)
        for i, targets in enumerate(edges):
            for j in targets:
                next_values[j] += values[i]
        values = next_values
        sequence.append(sum(values))
    return sequence, states, edges


def published_terms():
    path = Path(__file__).with_name("b204613.txt")
    terms = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        n, value = line.split()
        terms.append((int(n), int(value)))
    return terms


def main():
    sequence, states, edges = count_sequence(210)
    assert len(states) == 134
    assert sum(map(len, edges)) == 1114
    assert sequence[:5] == [907, 14199, 205117, 3084029, 45645661]

    terms = published_terms()
    assert [n for n, _ in terms] == list(range(1, 211))
    assert sequence == [value for _, value in terms]

    # OEIS uses a(n) = sum C[j-1] a(n-j), for n > 26.
    residuals = []
    for n in range(27, 211):
        residuals.append(
            sequence[n - 1]
            - sum(COEFFICIENTS[j - 1] * sequence[n - j - 1] for j in range(1, 24))
        )
    assert residuals == [0] * len(residuals)

    print({
        "states": len(states),
        "transitions": sum(map(len, edges)),
        "published_terms_checked": len(terms),
        "recurrence_interval": [27, 210],
        "zero_residuals": len(residuals),
    })
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
