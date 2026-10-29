"""Exact checks of the literal OEIS A185767 array predicate and certificate.

Run with Python 3. No recurrence is used to generate counts. The source data
are in b185767.txt; README.md gives their provenance and license.
"""

from collections import defaultdict
from itertools import product
import json
from pathlib import Path


COEFFICIENTS = (
    1, 98, 145, -1606, -1043, 9813, -11042, -14222, 67346,
    -40216, -64047, 36082, 29210, 8981, -23650, -12520, 7365, 448, -196,
)
SOURCE_PREFIX = (
    4456, 32456, 284256, 2680000, 27172072, 278758936, 2886249088,
    29857704328, 309306410040, 3201277096848, 33150509602128,
    343172969266424, 3553243907663632, 36786039676544160, 380867398106326120,
)
SOURCE_EXAMPLES = (
    ("10101110", "11010101", "00111011"),
    ("11101011", "11011101", "10111011"),
    ("11000100", "10101011", "01000111"),
)


def literal_valid(array):
    """Open rectangular boundary; ordinary integer main-diagonal trace."""
    height, width = len(array), len(array[0])
    traces = [
        [array[i][j] + array[i + 1][j + 1] for j in range(width - 1)]
        for i in range(height - 1)
    ]
    return (
        all(t[j] != t[j + 1] for t in traces for j in range(width - 2))
        and all(
            traces[i][j] != traces[i + 1][j]
            for i in range(height - 2)
            for j in range(width - 1)
        )
    )


def transfer(width):
    rows = list(product(range(2), repeat=width))
    states = []
    signatures = []
    by_first = defaultdict(list)
    for u in rows:
        for v in rows:
            t = tuple(u[j] + v[j + 1] for j in range(width - 1))
            if all(x != y for x, y in zip(t, t[1:])):
                by_first[u].append(len(states))
                states.append((u, v))
                signatures.append(t)
    edges = [
        [
            k for k in by_first[v]
            if all(x != y for x, y in zip(signatures[i], signatures[k]))
        ]
        for i, (u, v) in enumerate(states)
    ]
    return states, edges


def backward(edges, vector):
    return [sum(vector[j] for j in outgoing) for outgoing in edges]


def forward(edges, vector):
    result = [0] * len(edges)
    for i, outgoing in enumerate(edges):
        for j in outgoing:
            result[j] += vector[i]
    return result


def main():
    source = {}
    for line in Path(__file__).with_name("b185767.txt").read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            n, value = map(int, line.split())
            assert n not in source
            source[n] = value
    assert set(source) == set(range(1, 201))
    assert tuple(source[n] for n in range(1, 16)) == SOURCE_PREFIX

    states, edges = transfer(8)
    assert len(states) == 4456
    assert sum(map(len, edges)) == 32456

    # Independent literal-cell enumeration of all 2-by-8 arrays.
    brute = 0
    for bits in product(range(2), repeat=16):
        brute += literal_valid((bits[:8], bits[8:]))
    assert brute == source[1] == len(states)
    for example in SOURCE_EXAMPLES:
        assert literal_valid(tuple(tuple(map(int, row)) for row in example))

    # Count the first four source objects in the transposed orientation.
    transposed = []
    for n in range(1, 5):
        small_states, small_edges = transfer(n + 1)
        y = [1] * len(small_states)
        for _ in range(6):  # Eight rows, hence seven states and six edges.
            y = forward(small_edges, y)
        transposed.append(sum(y))
        assert transposed[-1] == source[n]

    powers = [[1] * len(states)]
    for _ in range(26):
        powers.append(backward(edges, powers[-1]))
    residual_vectors = []
    residuals = []
    for s in range(8):
        z = [
            powers[s + 19][i] - sum(
                c * powers[s + 18 - j][i]
                for j, c in enumerate(COEFFICIENTS)
            )
            for i in range(len(states))
        ]
        residual_vectors.append(z)
        residuals.append({
            "s": s, "nonzero": sum(x != 0 for x in z),
            "sum": sum(z), "max_abs": max(map(abs, z)),
        })
    assert sum(residual_vectors[5]) == sum(residual_vectors[6]) == 0
    assert all(x == 0 for x in residual_vectors[7])

    # Horner evaluation does not use the stored direct powers.
    horner = [1] * len(states)
    for c in COEFFICIENTS:
        horner = [x - c for x in backward(edges, horner)]
    for direct in residual_vectors:
        assert horner == direct
        horner = backward(edges, horner)

    # Independent forward path multiplicities through all 200 source terms.
    y = [1] * len(states)
    counts = {}
    for n in range(1, 201):
        counts[n] = sum(y)
        assert counts[n] == source[n], (n, counts[n], source[n])
        if n <= 27:
            assert counts[n] == sum(powers[n - 1])
        y = forward(edges, y)

    def scalar_residual(n):
        return counts[n] - sum(
            c * counts[n - j] for j, c in enumerate(COEFFICIENTS, start=1)
        )

    assert scalar_residual(24) == 78336000
    assert all(scalar_residual(n) == 0 for n in range(25, 201))
    output = {
        "states": len(states), "edges": sum(map(len, edges)),
        "coefficients": COEFFICIENTS, "residuals": residuals,
        "initial_counts_1_to_27": [counts[n] for n in range(1, 28)],
        "literal_2_by_8_count": brute, "source_examples_valid": 3,
        "transposed_counts_1_to_4": transposed,
        "source_bfile_comparison": "200/200 exact",
        "horner_comparison": "8/8 vectors exact",
        "scalar_recurrence_checks": "176/176 exact (25 <= n <= 200)",
        "source_threshold_residuals": {
            str(n): scalar_residual(n) for n in (24, 25, 26)
        },
    }
    print(json.dumps(output, indent=2))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
