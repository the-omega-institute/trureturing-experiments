#!/usr/bin/env python3
"""Exact finite checks for the literal OEIS A204614 array recurrence.

Run with Python 3; no third-party packages or network access are required.
The unbounded Cayley--Hamilton argument is explained in README.md.
"""

from collections import Counter
from itertools import product
from pathlib import Path
import hashlib
import json
import re


HERE = Path(__file__).resolve().parent


def determinant(u, v):
    return tuple(u[j] * v[j + 1] - u[j + 1] * v[j] for j in range(6))


def horizontal(d):
    return all(d[j] <= d[j + 1] for j in range(5))


def main():
    data = json.loads((HERE / "source-data.json").read_text())
    formula = next(line for line in data["entry_text"].splitlines()
                   if line.startswith("%F A204614 "))
    matches = re.findall(r"([+-]?\s*\d+)\*a\(n-(\d+)\)", formula)
    assert len(matches) == 43
    assert [int(i) for _, i in matches] == list(range(1, 44))
    coefficients = [int(c.replace(" ", "")) for c, _ in matches]
    assert formula.endswith("for n>47")
    assert "%O A204614 1,1" in data["entry_text"]

    rows = list(product((0, 1), repeat=7))
    pairs = {(u, v): determinant(u, v) for u in rows for v in rows
             if horizontal(determinant(u, v))}
    multiplicities = Counter((v, d) for (u, v), d in pairs.items())
    states = sorted(multiplicities)
    index = {s: i for i, s in enumerate(states)}
    initial = [multiplicities[s] for s in states]
    edges = []
    for v, d in states:
        edges.append(sorted(index[w, e] for w in rows
                            if (v, w) in pairs
                            for e in [pairs[v, w]]
                            if all(x <= y for x, y in zip(d, e))))
    dimension = len(states)
    assert dimension == 265
    assert sum(initial) == len(pairs) == 2470
    assert sum(map(len, edges)) == 3051
    print("TRANSFER_OK rows=128 states=265 initial_pairs=2470 edges=3051")

    # Independent integer-bit representation, including every potential edge.
    def bit_minor(u, v, j):
        return ((u >> j) & 1) * ((v >> (j + 1)) & 1) - (
            (u >> (j + 1)) & 1) * ((v >> j) & 1)

    def code(v):
        return sum(x << j for j, x in enumerate(v))

    bit_pairs = {}
    for u in range(128):
        for v in range(128):
            d = tuple(bit_minor(u, v, j) for j in range(6))
            if all(d[j] <= d[j + 1] for j in range(5)):
                bit_pairs[u, v] = d
    assert {(code(u), code(v)): d for (u, v), d in pairs.items()} == bit_pairs
    for i, (v, d) in enumerate(states):
        assert initial[i] == sum(bit_pairs.get((u, code(v))) == d
                                 for u in range(128))
        expected = [j for j, (w, e) in enumerate(states)
                    if bit_pairs.get((code(v), code(w))) == e
                    and all(d[k] <= e[k] for k in range(6))]
        assert edges[i] == expected
    three_rows = sum(all(d[j] <= e[j] for j in range(6))
                     for (u, v), d in bit_pairs.items()
                     for w in range(128)
                     if (e := bit_pairs.get((v, w))) is not None)
    assert three_rows == 62031
    print("INDEPENDENT_BIT_CHECK_OK all_states_weights_edges; three_rows=62031")

    def prefix_step(v):
        result = [0] * dimension
        for i, js in enumerate(edges):
            for j in js:
                result[j] += v[i]
        return result

    def suffix_step(v):
        return [sum(v[j] for j in js) for js in edges]

    # Prefix vectors count actual row extensions, without using the recurrence.
    counts = []
    prefix = initial[:]
    for _ in range(312):
        counts.append(sum(prefix))
        prefix = prefix_step(prefix)
    assert counts[:len(data["displayed_terms"])] == data["displayed_terms"]
    assert [n for n, _ in data["bfile_terms"]] == list(range(1, 211))
    for n, value in data["bfile_terms"]:
        assert counts[n - 1] == value, ("source fidelity", n)
    assert counts[1] == three_rows
    print("SOURCE_FIDELITY_OK displayed_terms=14 bfile_terms=210 indices=1..210")

    # T has source rows, target columns; b = T^4 P(T) 1.
    suffixes = [[1] * dimension]
    for _ in range(47):
        suffixes.append(suffix_step(suffixes[-1]))
    for n in range(1, 49):
        assert sum(x * y for x, y in zip(initial, suffixes[n - 1])) == counts[n - 1]
    b = [suffixes[47][j] - sum(coefficients[i - 1] * suffixes[47 - i][j]
                              for i in range(1, 44))
         for j in range(dimension)]
    prefix = initial[:]
    residuals = []
    for k in range(dimension):
        n = k + 48
        residual = counts[n - 1] - sum(coefficients[i - 1] * counts[n - i - 1]
                                      for i in range(1, 44))
        assert sum(x * y for x, y in zip(prefix, b)) == residual
        assert residual == 0, ("recurrence", n, residual)
        residuals.append(residual)
        prefix = prefix_step(prefix)
    boundary = counts[46] - sum(coefficients[i - 1] * counts[46 - i]
                                for i in range(1, 44))
    assert boundary == -79321900236876546048
    print("EXACT_RESIDUAL_BASE_OK 265_zero_residuals indices=48..312")
    print("PREFIX_SUFFIX_AGREEMENT_OK b=T^4*P(T)*1; r(k)=alpha*T^k*b")
    print("STRICT_BOUNDARY_OK residual_at_n47=-79321900236876546048")

    certificate = {
        "width": 7, "coefficients": coefficients, "states": states,
        "initial": initial, "edges": edges, "residual_vector_b": b,
        "counts_n1_to_n312": counts, "recurrence_residuals_n48_to_n312": residuals,
        "bfile_terms_matched": 210,
    }
    content = json.dumps(certificate, separators=(",", ":")) + "\n"
    (HERE / "certificate.json").write_text(content)
    print("CERTIFICATE_SHA256 " + hashlib.sha256(content.encode()).hexdigest())
    print("SCOPE exact_finite_checks; unbounded argument and Lean proof not executed")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
