#!/usr/bin/env python3
"""Exact finite checks for arXiv:2609.26127v1, Remarks 3.8 and 5.5.

Sparse integer polynomials and the source's Macaulay determinant quotient.
No floating point, numerical root finding, or assumed triangular resultant.
"""
from functools import lru_cache
from itertools import combinations, product
from math import factorial, prod


def add(p, q):
    r = dict(p)
    for a, c in q.items():
        r[a] = r.get(a, 0) + c
        if not r[a]:
            del r[a]
    return r


def scale(p, c):
    return {a: c * v for a, v in p.items() if c * v}


def mul(p, q):
    r = {}
    for a, c in p.items():
        for b, v in q.items():
            z = tuple(x + y for x, y in zip(a, b))
            r[z] = r.get(z, 0) + c * v
    return {a: c for a, c in r.items() if c}


def xp(n, i, d):
    return {tuple(d if j == i else 0 for j in range(n)): 1}


@lru_cache(None)
def monomials(n, d):
    if n == 0:
        return ((),) if d == 0 else ()
    return tuple((i,) + a for i in range(d + 1)
                 for a in monomials(n - 1, d - i))


def complete(n, indices, degree):
    r = {}
    for a in monomials(len(indices), degree):
        b = [0] * n
        for i, e in zip(indices, a):
            b[i] += e
        r[tuple(b)] = 1
    return r


def divide_difference(p, i, j):
    """Exact division by X_i-X_j, using synthetic polynomial division."""
    r, q = dict(p), {}
    while r:
        a = max(r, key=lambda v: v[i])
        assert a[i] > 0, ("nonzero remainder", r)
        b = list(a)
        b[i] -= 1
        b = tuple(b)
        c = r[a]
        q = add(q, {b: c})
        other = list(b)
        other[j] += 1
        r = add(r, {a: -c})
        r = add(r, {tuple(other): c})
    return q


def divided(F, indices):
    if len(indices) == 1:
        return F[indices[0]]
    p = add(divided(F, indices[:-1]),
            scale(divided(F, indices[:-2] + indices[-1:]), -1))
    return divide_difference(p, indices[-2], indices[-1])


def specialize(p, mapping, n):
    r = {}
    for a, c in p.items():
        b = [0] * n
        for i, e in enumerate(a):
            b[mapping[i]] += e
        r = add(r, {tuple(b): c})
    return r


def derivative(p, i):
    r = {}
    for a, c in p.items():
        if a[i]:
            b = list(a)
            b[i] -= 1
            r[tuple(b)] = c * a[i]
    return r


def bareiss(a):
    """Fraction-free integer determinant, with exact divisibility asserted."""
    a = [row[:] for row in a]
    n, sign, previous = len(a), 1, 1
    if not n:
        return 1
    for k in range(n - 1):
        pivot = next((i for i in range(k, n) if a[i][k]), None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        v = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                t = v * a[i][j] - a[i][k] * a[k][j]
                assert t % previous == 0
                a[i][j] = t // previous
            a[i][k] = 0
        previous = v
    return sign * a[-1][-1]


MAX_MATRIX = 0
DETERMINANTS = 0
DENSE_CROSSCHECKS = 0


def determinant(entries, n):
    global DETERMINANTS, DENSE_CROSSCHECKS
    DETERMINANTS += 1
    triangular = (all(i <= j for i, j in entries)
                  or all(i >= j for i, j in entries))
    value = prod(entries.get((i, i), 0) for i in range(n)) if triangular else None
    if not triangular or n <= 60:
        dense = [[entries.get((i, j), 0) for j in range(n)] for i in range(n)]
        v = bareiss(dense)
        if triangular:
            assert value == v
            DENSE_CROSSCHECKS += 1
        value = v
    return value


def resultant(F, degrees):
    """Source §2.2: det(M)/det(D), with the same row/column ordering."""
    global MAX_MATRIX
    n = len(F)
    assert n > 0 and all(e > 0 for e in degrees)
    for p, e in zip(F, degrees):
        assert all(len(a) == n and sum(a) == e for a in p)
    t = sum(e - 1 for e in degrees) + 1
    mons = monomials(n, t)
    index = {a: i for i, a in enumerate(mons)}
    MAX_MATRIX = max(MAX_MATRIX, len(mons))
    entries, dod = {}, []
    for col, beta in enumerate(mons):
        divisors = [i for i, e in enumerate(degrees) if beta[i] >= e]
        i = min(divisors)
        if len(divisors) >= 2:
            dod.append(col)
        shift = list(beta)
        shift[i] -= degrees[i]
        for a, c in F[i].items():
            row = index[tuple(x + y for x, y in zip(a, shift))]
            entries[row, col] = c
    subindex = {v: i for i, v in enumerate(dod)}
    sub = {(subindex[i], subindex[j]): c for (i, j), c in entries.items()
           if i in subindex and j in subindex}
    numerator = determinant(entries, len(mons))
    denominator = determinant(sub, len(dod))
    assert denominator != 0, "Macaulay quotient specialization is indeterminate"
    assert numerator % denominator == 0
    return numerator // denominator


@lru_cache(None)
def partitions(n, maximum=None):
    maximum = n if maximum is None else min(n, maximum)
    if n == 0:
        return ((),)
    return tuple((i,) + p for i in range(maximum, 0, -1)
                 for p in partitions(n - i, i))


def multiplicity(p):
    denominator = prod(factorial(i) for i in p)
    denominator *= prod(factorial(p.count(i)) for i in set(p))
    assert factorial(sum(p)) % denominator == 0
    return factorial(sum(p)) // denominator


def falling(d, k):
    return prod(range(d - k + 1, d + 1))


def m0(n, d):
    total = 0
    for p in partitions(n):
        k = len(p)
        if k <= d:
            total += multiplicity(p) * sum(falling(d, k) // (d - j)
                                           for j in range(k))
    value = n * d ** (n - 1) - total
    assert value >= 0
    assert (value == 0) == (n <= d)
    return value


def compositions(n):
    if n == 0:
        return ((),)
    return tuple((i,) + p for i in range(1, n + 1)
                 for p in compositions(n - i))


def fidelity():
    # Source §4, displayed d=3,n=4,p=q=2 example: first binary resultant.
    x = [xp(4, i, 1) for i in range(4)]
    e1 = complete(4, tuple(range(4)), 1)
    e2, e3 = {}, {}
    for size, target in [(2, e2), (3, e3)]:
        for inds in combinations(range(4), size):
            term = {tuple(0 for _ in range(4)): 1}
            for i in inds:
                term = mul(term, x[i])
            target.update(add(target, term))
    cases = [(1, 0, 0), (1, 1, 1), (3, 2, 1), (2, -1, 2), (-2, 3, -1)]
    values = []
    for a, b, c in cases:
        p = add(add(scale(e3, a), scale(mul(e2, e1), b)),
                scale(mul(mul(e1, e1), e1), c))
        F = [specialize(derivative(p, i), [0, 0, 1, 1], 2) for i in [0, 2]]
        actual = resultant(F, [2, 2])
        expected = -3 * (2 * b + a) ** 3 * (16 * c + 6 * b + a)
        assert actual == expected, (a, b, c, actual, expected)
        values.append(actual)
    print("FIDELITY source_section=4 degree=3 n=4 p=q=2 binary_factor_cases=5 values=" + str(values))


def small_cases():
    assert resultant([xp(2, 0, 1), xp(2, 1, 1)], [1, 1]) == 1
    assert resultant([xp(2, 1, 1), xp(2, 0, 1)], [1, 1]) == -1
    tuples = factors = extracted = dd_checks = 0
    shapes = set()
    cases = [(n, sizes, d) for n in range(2, 6)
             for sizes in compositions(n) for d in range(1, 4)]
    cases += [(6, (2, 2, 2), d) for d in [1, 2]]
    for n, sizes, d in cases:
        tuples += 1
        F = [xp(n, i, d) for i in range(n)]
        lhs = resultant(F, [d] * n)
        assert lhs == 1
        prefactor = 1
        for s, size in enumerate(sizes):
            start = sum(sizes[:s])
            if size > d:
                constant = divided(F, tuple(range(start, start + d + 1)))
                assert constant == {(0,) * n: 1}
                mu = d ** (n - size) * m0(size, d)
                prefactor *= constant[(0,) * n] ** mu
                extracted += 1
            else:
                assert m0(size, d) == 0
        rhs = prefactor
        allowed = [[p for p in partitions(size) if len(p) <= d] for size in sizes]
        for parts in product(*allowed):
            lengths = tuple(map(len, parts))
            N = sum(lengths)
            mapping, representatives = [], []
            orig_start, new_start = 0, 0
            for part in parts:
                reps, local = [], 0
                for j, count in enumerate(part):
                    reps.append(orig_start + local)
                    mapping.extend([new_start + j] * count)
                    local += count
                representatives.append(reps)
                orig_start += sum(part)
                new_start += len(part)
            chain, degrees = [], []
            offset = 0
            for reps in representatives:
                for j in range(1, len(reps) + 1):
                    q = specialize(divided(F, tuple(reps[:j])), mapping, N)
                    h = complete(N, tuple(range(offset, offset + j)), d - j + 1)
                    assert q == h
                    chain.append(q)
                    degrees.append(d - j + 1)
                    dd_checks += 1
                offset += len(reps)
            res = resultant(chain, degrees)
            assert res == 1
            shapes.add((lengths, d))
            rhs *= res ** prod(multiplicity(p) for p in parts)
            # Independently check the blockwise coefficient-degree identity (star_r).
            factors += 1
        assert lhs == rhs == 1
        for s, size in enumerate(sizes):
            degree_sum = 0
            for parts in product(*allowed):
                k = tuple(map(len, parts))
                D = prod(falling(d, j) for j in k)
                degree_sum += prod(multiplicity(p) for p in parts) * sum(
                    D // (d - j) for j in range(k[s]))
            mu = d ** (n - size) * m0(size, d)
            assert size * d ** (n - 1) == mu + degree_sum
    print(f"SMALL_CASES all_ordered_positive_block_compositions total_n=2..5 d=1..3 plus_(2,2,2)_d=1,2 tuples={tuples}")
    print(f"FACTORS partition_tuples={factors} distinct_chain_degree_shapes={len(shapes)} extracted_constants={extracted} divided_difference_checks={dd_checks}")
    print(f"MACAULAY determinants={DETERMINANTS} dense_Bareiss_crosschecks={DENSE_CROSSCHECKS} max_matrix_dimension={MAX_MATRIX}")
    print("SIGN_CONTROL simultaneous_identity=+1 swapped_linear_polynomials=-1")
    print("SCOPE finite_exact_checks_only; no assertion of a universal proof or completed Lean settlement")


def main():
    fidelity()
    small_cases()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
