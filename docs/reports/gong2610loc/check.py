#!/usr/bin/env python3
"""Exact finite checks for Gong et al., arXiv:2610.04666v1, Section V(c).

Independent implementation from the paper's formulas. This is numerical evidence
on finite cubes, not a proof of the unbounded preregistered statement.
"""
from fractions import Fraction
from itertools import combinations
from math import comb


def dependencies(values, ell):
    return sum(1 << j for j in range(ell)
               if any(values[s] != values[s ^ (1 << j)] for s in range(1 << ell)))


def bound(d, k):
    return 4 ** max(d - 2, 0) * sum(comb(k, j) for j in range(1, min(k, d) + 1))


def sign_table(table):
    return tuple(1 - 2 * b for b in table)


def walsh(ell, mask):
    return tuple(1 - 2 * ((s & mask).bit_count() % 2) for s in range(1 << ell))


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def all_families(ell):
    """All unordered balanced orthogonal families, including every output polarity.

    Distinct outputs cannot be equal or complements under orthogonality, so using
    increasing truth-table indices loses only output order, which changes no test.
    """
    size = 1 << ell
    tables = [tuple((mask >> s) & 1 for s in range(size))
              for mask in range(1 << size) if mask.bit_count() * 2 == size]
    signs = [sign_table(t) for t in tables]
    adj = [{j for j in range(i + 1, len(tables)) if dot(signs[i], signs[j]) == 0}
           for i in range(len(tables))]

    def visit(chosen, pool):
        yield tuple(tables[i] for i in chosen)
        for offset, i in enumerate(pool):
            yield from visit(chosen + [i], [j for j in pool[offset + 1:] if j in adj[i]])

    yield from visit([], list(range(len(tables))))


def exhaustive():
    totals = {'families': 0, 'pinned_families': 0, 'parameter_cases': 0,
              'eq6_outputs': 0, 'bessel_cases': 0, 'granularity_functions': 0}
    maxima = {}
    for ell in range(4):
        size = 1 << ell
        characters = [walsh(ell, a) for a in range(size)]
        # Orthogonality/completeness on each cube, as an independent character check.
        for a in range(size):
            for b in range(size):
                assert dot(characters[a], characters[b]) == (size if a == b else 0)
        local_counts = {'families': 0, 'pinned_families': 0, 'parameter_cases': 0}
        for family in all_families(ell):
            totals['families'] += 1
            local_counts['families'] += 1
            n = len(family)
            xs = [sign_table(t) for t in family]
            deps = [dependencies(t, ell) for t in family]
            locality = max((t.bit_count() for t in deps), default=0)
            weight = tuple(sum(t[s] for t in family) for s in range(size))
            coeff = [[Fraction(dot(x, c), size) for c in characters] for x in xs]
            for i, t in enumerate(deps):
                assert coeff[i][0] == 0
                assert sum(c * c for c in coeff[i]) == 1
                delta = Fraction(1, 2 ** max(t.bit_count() - 2, 0))
                for a, c in enumerate(coeff[i]):
                    if a & ~t:
                        assert c == 0
                    if c:
                        assert abs(c) >= delta
                totals['granularity_functions'] += 1
            for a in range(size):
                energy = sum((row[a] ** 2 for row in coeff), Fraction(0))
                residual = [Fraction(characters[a][s]) -
                            sum((coeff[i][a] * xs[i][s] for i in range(n)), Fraction(0))
                            for s in range(size)]
                assert sum(q * q for q in residual) / size == 1 - energy
                assert energy <= 1
                totals['bessel_cases'] += 1
            # n <= 2^ell-1 <= 7. Prime 11 represents every odd p > 7:
            # W is between 0 and n, so residue equality is then integer equality.
            for p in (3, 5, 7, 11):
                residue = tuple(w % p for w in weight)
                residue_dep = dependencies(residue, ell)
                for k_mask in range(size):
                    if residue_dep & ~k_mask:
                        continue
                    totals['pinned_families'] += 1
                    local_counts['pinned_families'] += 1
                    k = k_mask.bit_count()
                    for i, x in enumerate(xs):
                        # Source Eq. (6), evaluated both by fibers and pointwise.
                        fibers = {}
                        for s in range(size):
                            key = s & k_mask
                            fibers[key] = fibers.get(key, 0) + x[s]
                        grouped = sum(residue[sigma] * a for sigma, a in fibers.items())
                        assert grouped == sum(residue[s] * x[s] for s in range(size))
                        assert ell >= 1
                        assert (grouped + 2 ** (ell - 1)) % p == 0
                        assert any(a != 0 for a in fibers.values())
                        assert sum(weight[s] * x[s] for s in range(size)) == -size // 2
                        totals['eq6_outputs'] += 1
                    for d in range(locality, 4):
                        totals['parameter_cases'] += 1
                        local_counts['parameter_cases'] += 1
                        b = bound(d, k)
                        assert n <= b <= 4 ** max(d - 2, 0) * (2 ** k - 1)
                        delta = Fraction(1, 2 ** max(d - 2, 0))
                        selected = [a for a in range(1, size)
                                    if not (a & ~k_mask) and a.bit_count() <= d]
                        assert len(selected) == sum(comb(k, j) for j in range(1, min(k, d) + 1))
                        for row in coeff:
                            assert sum((row[a] ** 2 for a in selected), Fraction(0)) >= delta ** 2
                        assert n * delta ** 2 <= len(selected)
                        maxima[(d, k)] = max(maxima.get((d, k), 0), n)
        print('EXHAUSTIVE', 'ell=' + str(ell), local_counts)
    print('TOTALS', totals)
    print('MAXIMUM_N', sorted(maxima.items()))


def sharpness():
    count = 0
    for k in range(7):
        for d in (1, 2):
            masks = [a for a in range(1, 1 << k) if a.bit_count() <= d]
            family = [tuple((s & a).bit_count() % 2 for s in range(1 << k)) for a in masks]
            xs = [sign_table(t) for t in family]
            n = len(family)
            assert n == bound(d, k)
            assert all(sum(x) == 0 for x in xs)
            for i, j in combinations(range(n), 2):
                assert dot(xs[i], xs[j]) == 0
            assert all(dependencies(t, k).bit_count() <= d for t in family)
            weight = [sum(t[s] for t in family) for s in range(1 << k)]
            for p in (3, 5, 7, 11):
                r = [w % p for w in weight]
                assert all(0 <= a < p for a in r)
                count += 1
    assert bound(3, 3) == 28
    print('SHARPNESS', 'k=0..6 d=1,2 primes=3,5,7,11 cases=' + str(count))
    print('PARAMETER_BOUND', 'd=k=3 N<=28')


def fidelity_example():
    # All three nonempty characters on K={0,1}, one unused seed, p=3.
    ell, p, k_mask = 3, 3, 3
    family = [tuple((s & a).bit_count() % 2 for s in range(8)) for a in (1, 2, 3)]
    weight = [sum(t[s] for t in family) for s in range(8)]
    readings = []
    for t in family:
        x = sign_table(t)
        fibers = [sum(x[s] for s in range(8) if s & k_mask == sigma) for sigma in range(4)]
        lhs = sum((weight[sigma] % p) * fibers[sigma] for sigma in range(4))
        assert lhs == -4
        assert lhs % p == (-2 ** (ell - 1)) % p
        readings.append({'fiber_sums': fibers, 'eq6_lhs': lhs, 'residue_mod_3': lhs % p})
    print('SOURCE_EQ6', readings)


if __name__ == '__main__':
    fidelity_example()
    exhaustive()
    sharpness()
    print('ALL_CHECKS_OK')
