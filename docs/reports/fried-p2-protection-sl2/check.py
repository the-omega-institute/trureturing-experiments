#!/usr/bin/env python3
"""Exact finite checks for Fried's Section 10 CA-group protection target.

The family argument belongs to the preregistration. No floating point is used.
The q=64 calculation does not construct its full commuting matrix.
"""

from collections import Counter, defaultdict
from itertools import product
import json

import sympy
from sympy import Matrix, ZZ, factorint
from sympy.matrices.normalforms import smith_normal_form


def poly_rem(a, b):
    while a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def binary_field(m, modulus):
    assert modulus.bit_length() == m + 1
    # Any reducible degree-m polynomial has a monic factor of degree <= m//2.
    for d in range(1, m // 2 + 1):
        for divisor in range(1 << d, 1 << (d + 1)):
            assert poly_rem(modulus, divisor) != 0
    q = 1 << m

    def mul(a, b):
        out = 0
        while b:
            if b & 1:
                out ^= a
            b >>= 1
            a <<= 1
            if a & q:
                a ^= modulus
        return out

    table = [[mul(a, b) for b in range(q)] for a in range(q)]
    inverse = [0] + [next(b for b in range(1, q) if table[a][b] == 1)
                     for a in range(1, q)]
    for a, b, c in product(range(q), repeat=3):
        assert table[table[a][b]][c] == table[a][table[b][c]]
        assert table[a][b ^ c] == table[a][b] ^ table[a][c]
    return q, table, inverse


def matmul(a, b, table, add):
    u, v, w, x = a
    y, z, t, s = b
    return (add(table[u][y], table[v][t]),
            add(table[u][z], table[v][s]),
            add(table[w][y], table[x][t]),
            add(table[w][z], table[x][s]))


def sl2_elements(q, table, inverse, add, neg):
    out = []
    for a, b in product(range(q), repeat=2):
        if a:
            for c in range(q):
                d = table[inverse[a]][add(1, table[b][c])]
                out.append((a, b, c, d))
        elif b:
            c = neg(inverse[b])
            out.extend((0, b, c, d) for d in range(q))
    assert len(out) == q * (q * q - 1)
    assert len(set(out)) == len(out)
    for a, b, c, d in out:
        assert add(table[a][d], neg(table[b][c])) == 1
    return out


def fidelity_sl2_3():
    """Reproduce the source's Table 1: SL(2,3), b=7, Smith torsion 6."""
    q = 3
    table = [[a * b % q for b in range(q)] for a in range(q)]
    inv = [0, 1, 2]
    add = lambda a, b: (a + b) % q
    neg = lambda a: -a % q
    group = sl2_elements(q, table, inv, add, neg)
    mul = lambda a, b: matmul(a, b, table, add)
    centralizers = [frozenset(j for j, h in enumerate(group) if mul(g, h) == mul(h, g))
                   for g in group]
    center = {i for i, c in enumerate(centralizers) if len(c) == len(group)}
    blocks = set(centralizers) - {frozenset(range(len(group)))}
    assert len(center) == 2 and len(blocks) == 7
    for block in blocks:
        for i, j in product(block, repeat=2):
            assert mul(group[i], group[j]) == mul(group[j], group[i])
    base = min(center)
    rows = [[int(j in centralizers[i]) - int(j in centralizers[base])
             for j in range(len(group))] for i in range(len(group)) if i != base]
    smith = smith_normal_form(Matrix(rows), domain=ZZ)
    factors = [abs(int(smith[i, i])) for i in range(min(smith.shape)) if smith[i, i]]
    assert factors == [1] * 6 + [6]
    print(json.dumps({"fidelity": "arXiv:2607.05698v2, Table 1, SL(2,3)",
                      "order": 24, "center_order": 2, "b": 7,
                      "nonzero_smith_factors": factors, "conductor": [2, 3]}, sort_keys=True))


def field_rank(rows, table, inverse):
    rows = [list(row) for row in rows]
    rank = 0
    for col in range(len(rows[0])):
        pivot = next((r for r in range(rank, len(rows)) if rows[r][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        scale = inverse[rows[rank][col]]
        rows[rank] = [table[scale][x] for x in rows[rank]]
        for r in range(len(rows)):
            if r != rank:
                scale = rows[r][col]
                rows[r] = [x ^ table[scale][y] for x, y in zip(rows[r], rows[rank])]
        rank += 1
    return rank


def check_even_sl2(m, modulus):
    q, table, inverse = binary_field(m, modulus)
    add = lambda a, b: a ^ b
    group = sl2_elements(q, table, inverse, add, lambda a: a)
    identity = (1, 0, 0, 1)
    mul = lambda a, b: matmul(a, b, table, add)
    U = lambda t: (1, t, 0, 1)
    L = lambda t: (1, 0, t, 1)

    def sig(g):
        a, b, c, d = g
        vector = (b, c, a ^ d)  # Coordinates in Mat_2(F)/FI.
        pivot = next((x for x in vector if x), None)
        return None if pivot is None else tuple(table[inverse[pivot]][x] for x in vector)

    blocks = defaultdict(list)
    for g in group:
        if g != identity:
            blocks[sig(g)].append(g)
    assert None not in blocks
    assert len(blocks) == q * q + q + 1
    expected_signatures = {(1, a, b) for a, b in product(range(q), repeat=2)}
    expected_signatures |= {(0, 1, a) for a in range(q)} | {(0, 0, 1)}
    assert set(blocks) == expected_signatures

    # Independently check that every representative has a 2-dimensional commutant.
    # For A=(a,b;c,d), coefficients of XA-AX in the four entries of X are:
    for block in blocks.values():
        a, b, c, d = representative = block[0]
        equations = [[0, c, b, 0], [b, a ^ d, 0, b],
                     [c, 0, a ^ d, c], [0, c, b, 0]]
        assert field_rank(equations, table, inverse) == 2
        for g in block:
            assert mul(representative, g) == mul(g, representative)

    pair_checks = 0
    if q <= 8:
        signatures = {g: sig(g) for g in group}
        for g, h in product(group, repeat=2):
            assert (mul(g, h) == mul(h, g)) == (
                g == identity or h == identity or signatures[g] == signatures[h])
            pair_checks += 1

    # Each transvection is a displayed diagonal commutator. Characteristic two
    # turns inverses of U,L into themselves and a^2-1 into a^2+1.
    a = 2
    coefficient = table[a][a] ^ 1
    assert coefficient != 0
    diagonal = (a, 0, 0, inverse[a])
    diagonal_inv = (inverse[a], 0, 0, a)
    lower_coefficient = table[inverse[a]][inverse[a]] ^ 1
    assert lower_coefficient != 0
    for t in range(q):
        s = table[inverse[coefficient]][t]
        assert mul(mul(mul(diagonal, U(s)), diagonal_inv), U(s)) == U(t)
        s = table[inverse[lower_coefficient]][t]
        assert mul(mul(mul(diagonal, L(s)), diagonal_inv), L(s)) == L(t)

    # Verify an actual U,L factorization for every matrix, including pivot zero.
    swap = mul(mul(U(1), L(1)), U(1))
    for original in group:
        g = original
        prefix = identity
        if g[0] == 0:
            prefix, g = swap, mul(swap, g)
        a, b, c, d = g
        factors = [L(table[c][inverse[a]]), U(a), L(inverse[a]), U(a),
                   U(1), L(1), U(1), U(table[b][inverse[a]])]
        reconstructed = prefix
        for factor in factors:
            reconstructed = mul(reconstructed, factor)
        assert reconstructed == original

    types = Counter(len(block) + 1 for block in blocks.values())
    assert types == Counter({q: q + 1, q - 1: q * (q + 1) // 2,
                             q + 1: q * (q - 1) // 2})
    b = len(blocks)
    assert b - 1 == q * (q + 1)
    record = {"q": q, "binary_modulus": modulus, "order": len(group), "b": b,
              "centralizer_orders_and_counts": dict(sorted(types.items())),
              "commutant_rank_checks": b, "commuting_pair_checks": pair_checks,
              "generation_factorizations": len(group),
              "b_minus_one_factorization": {str(k): int(v) for k, v in factorint(b - 1).items()}}
    if q == 64:
        assert len(group) % 9 == 0 and (b - 1) % 3 == 2
        record.update({"p": 3, "p_squared_divides_perfect_group_order": True,
                       "p_divides_b_minus_one": False,
                       "order_factorization": {str(k): int(v) for k, v in factorint(len(group)).items()}})
    print(json.dumps(record, sort_keys=True))


def family_arithmetic_samples():
    checked = []
    for p, k in product((3, 5, 7, 11, 13), (1, 2)):
        assert sympy.isprime(p)
        exponent = p * (p - 1) * k
        q = 1 << exponent
        b = q * q + q + 1
        order = q * (q * q - 1)
        assert q >= 64 and q % (p * p) == 1
        assert order % (p * p) == 0 and (b - 1) % p == 2
        checked.append({"p": p, "k": k, "exponent": exponent})
    print(json.dumps({"arithmetic_samples_only": checked}, sort_keys=True))


if __name__ == "__main__":
    print("Exact integer and finite-field checks; sympy=" + sympy.__version__)
    fidelity_sl2_3()
    for m, modulus in ((2, 0b111), (3, 0b1011), (6, 0b1000011)):
        check_even_sl2(m, modulus)
    family_arithmetic_samples()
    print("ALL_CHECKS_OK")
