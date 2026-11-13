#!/usr/bin/env python3
"""Exact, finite checks for Ma's size-four distinct-powers question.

Independent implementation; Python standard library only. No finite test here
establishes the universal finite-group statement.
"""
from collections import Counter
from itertools import combinations, permutations, product
from math import gcd

L0 = ((0, 1, 2, 3), (1, 2, 3, 0), (2, 3, 0, 1), (3, 0, 1, 2))
L1 = ((0, 1, 2, 3), (1, 2, 0, 4), (2, 0, 3, 1), (4, 3, 1, 0))


def transversal(array):
    n = len(array)
    return next((p for p in permutations(range(n))
                 if len({array[i][p[i]] for i in range(n)}) == n), None)


def normalize(flat):
    names = {}
    return tuple(names.setdefault(x, len(names)) for x in flat)


def pattern_orbit(array):
    n = len(array)
    out = set()
    for a in (array, tuple(zip(*array))):
        for r in permutations(range(n)):
            for c in permutations(range(n)):
                out.add(normalize(a[i][j] for i in r for j in c))
    return out


def latin_obstructions(n):
    """All Latin arrays without transversals, modulo symbol renaming only.

    Restricted-growth labelling gives first row 0,...,n-1. Enumerate the
    first n-1 rows. A last-row cell must belong to every rainbow partial
    transversal's symbol set on the complementary columns. This is necessary
    and sufficient to block every full transversal; complete all allowed
    last rows with the Latin and restricted-growth conditions.
    """
    first = tuple(range(n))
    rows = [first]
    column_masks = [1 << x for x in first]
    obstructions = set()
    prefixes = 0
    assignment_nodes = 0

    def finish(maximum):
        nonlocal prefixes
        prefixes += 1
        allowed = []
        for missing in range(n):
            cols = tuple(j for j in range(n) if j != missing)
            intersection = None
            for p in permutations(cols):
                symbols = {rows[i][p[i]] for i in range(n - 1)}
                if len(symbols) == n - 1:
                    intersection = symbols if intersection is None else intersection & symbols
            # The finite n*n bound also covers a fresh last-row symbol.
            allowed.append(set(range(n*n)) if intersection is None else intersection)
        if any(not a for a in allowed):
            return

        def last(col, row, used, high):
            if col == n:
                array = tuple(rows) + (tuple(row),)
                assert transversal(array) is None
                obstructions.add(tuple(x for r in array for x in r))
                return
            for v in sorted(allowed[col]):
                if v > high + 1 or (used | column_masks[col]) & (1 << v):
                    continue
                last(col + 1, row + [v], used | (1 << v), max(high, v))
        last(0, [], 0, maximum)

    def extend(col, row, used, high):
        nonlocal assignment_nodes
        assignment_nodes += 1
        if col == n:
            rows.append(tuple(row))
            for j, v in enumerate(row):
                column_masks[j] |= 1 << v
            if len(rows) == n - 1:
                finish(high)
            else:
                extend(0, [], 0, high)
            for j, v in enumerate(row):
                column_masks[j] ^= 1 << v
            rows.pop()
            return
        for v in range(high + 2):
            if (used | column_masks[col]) & (1 << v):
                continue
            extend(col + 1, row + [v], used | (1 << v), max(high, v))

    if n == 2:
        finish(n - 1)
    else:
        extend(0, [], 0, n - 1)
    return obstructions, prefixes, assignment_nodes


def classification_checks():
    obs3, p3, nodes3 = latin_obstructions(3)
    assert not obs3
    print(f'LATIN order=3 prefixes={p3} assignment_nodes={nodes3} transversal_free=0')
    obs4, p4, nodes4 = latin_obstructions(4)
    orbit0, orbit1 = pattern_orbit(L0), pattern_orbit(L1)
    assert transversal(L0) is None and transversal(L1) is None
    assert not orbit0 & orbit1
    assert obs4 == orbit0 | orbit1
    print(f'LATIN order=4 prefixes={p4} assignment_nodes={nodes4} '
          f'transversal_free={len(obs4)} L0={len(orbit0)} L1={len(orbit1)} classes=2')


def faithful_checks():
    # Ma, Remark 2.2: a multiset, not four distinct units.
    values = (1, 1, 1, 10)
    bad = sum(len({a * u % 11 for a, u in zip((4, 5, 6, 7), p)}) < 4
              for p in permutations(values))
    assert bad == 24
    print('FIDELITY Remark_2.2 F11 block=4,5,6,7 multiset=1,1,1,-1 bad_labelled_assignments=24/24')
    # Ma, Lemma 4.1: the six bad four-entry multisets over F7.
    # Enumerate sorted multisets by combinations with repetition, independently.
    from itertools import combinations_with_replacement
    bad7 = {u for u in combinations_with_replacement(range(1, 7), 4)
            if transversal(tuple(tuple(a * v % 7 for v in u) for a in range(1, 5))) is None}
    expected = {tuple(sorted((s, s, 2*s % 7, 2*s % 7))) for s in range(1, 7)}
    assert bad7 == expected and len(bad7) == 6
    print('FIDELITY Lemma_4.1 F7 bad_multisets=' + repr(sorted(bad7)))


def local_checks():
    # Distinct units modulo full m, not necessarily after prime reduction.
    # Includes prime squares (especially 49 and 121) with repeated reductions.
    count = 0
    repeated_reductions = 0
    for m in tuple(range(2, 32)) + (49, 121):
        prime_divs = [p for p in range(2, m+1) if m % p == 0 and
                      all(p % d for d in range(2, p))]
        pmin = min(prime_divs)
        if pmin <= 5:
            continue
        units = tuple(u for u in range(1, m) if gcd(u, m) == 1)
        if m == 121:
            units = (1, 2, 3, 4, 12, 23, 34, 45)
        for u in combinations(units, 4):
            if len({v % pmin for v in u}) < 4:
                repeated_reductions += 1
            for r in range(1, min(12, pmin - 4)):
                # n=r+3 satisfies p > n+1 exactly when p > r+4.
                array = tuple(tuple(a * v % m for v in u) for a in range(r, r+4))
                assert transversal(array) is not None, (m, r, u)
                count += 1
    print(f'FOUR_BLOCK m<=31_or_m=49_all_units_or_m=121_eight_named_units r<=11 checked={count} unit_sets_with_repeated_prime_reductions={repeated_reductions}')
    assert repeated_reductions > 0
    # Exact power sums and substitution in characteristic not 2.
    from fractions import Fraction
    r = Fraction(-3, 2)
    assert 4*r + 6 == 0 and 4*r*r + 12*r + 14 == 5
    print('POWER_SUMS at_r=-3/2 first=0 second=5 exact_rational=True')


def group_data(elements, multiply, identity):
    index = {x: i for i, x in enumerate(elements)}
    table = [[index[multiply(x, y)] for y in elements] for x in elements]
    one = index[identity]
    # Exact associativity check for these small test groups.
    for x in range(len(elements)):
        assert table[one][x] == x == table[x][one]
        assert any(table[x][y] == one == table[y][x] for y in range(len(elements)))
        for y in range(len(elements)):
            for z in range(len(elements)):
                assert table[table[x][y]][z] == table[x][table[y][z]]
    powers = []
    orders = []
    subgroups = []
    for x in range(len(elements)):
        row, y = [], one
        for _ in range(9):
            y = table[y][x]
            row.append(y)
        powers.append(row)
        cyclic = {one}
        y = x
        while y != one:
            assert y not in cyclic
            cyclic.add(y)
            y = table[y][x]
        orders.append(len(cyclic))
        subgroups.append(frozenset(cyclic))
    return powers, orders, subgroups


def has_assignment(subset, powers):
    # Complete backtracking; a returned solution is checked independently below.
    n = len(subset)
    assignment = [-1] * n
    def visit(i, remaining, used):
        if i == n:
            return tuple(assignment)
        for x in remaining:
            y = powers[x][i]
            if y not in used:
                assignment[i] = x
                result = visit(i+1, tuple(z for z in remaining if z != x), used | {y})
                if result is not None:
                    return result
        return None
    result = visit(0, subset, set())
    if result is not None:
        assert sorted(result) == sorted(subset)
        assert len({powers[x][i] for i, x in enumerate(result)}) == n
    return result is not None


def group_checks():
    groups = []
    for m in range(1, 32, 2):
        groups.append((f'C{m}', tuple(range(m)), lambda x, y, m=m: (x+y) % m, 0, None))
    for p in (3, 5):
        elems = tuple(product(range(p), repeat=2))
        groups.append((f'C{p}xC{p}', elems, lambda x, y, p=p: tuple((a+b) % p for a,b in zip(x,y)), (0,0), None))
    # C7 semidirect C3, generator of C3 acts by multiplication by 2 modulo 7.
    elems = tuple(product(range(7), range(3)))
    groups.append(('C7_semidirect_C3', elems,
                   lambda x,y: ((x[0] + pow(2,x[1],7)*y[0]) % 7, (x[1]+y[1]) % 3), (0,0), None))
    elems = tuple(product(range(9), range(3)))
    groups.append(('C9xC3', elems, lambda x,y: ((x[0]+y[0])%9, (x[1]+y[1])%3), (0,0), None))
    groups.append(('C9_semidirect_C3', elems,
                   lambda x,y: ((x[0] + pow(4,x[1],9)*y[0])%9, (x[1]+y[1])%3), (0,0), None))
    elems = tuple(product(range(3), repeat=3))
    groups.append(('C3xC3xC3', elems, lambda x,y: tuple((a+b)%3 for a,b in zip(x,y)), (0,0,0), None))
    groups.append(('Heisenberg3', elems, lambda x,y: ((x[0]+y[0])%3, (x[1]+y[1])%3,
                   (x[2]+y[2]+x[0]*y[1])%3), (0,0,0), None))
    # This explicitly bounded 9-point domain supplies nonvacuous n=6,...,9
    # size-four-class cases in the actual 121-element group.
    elems = tuple(product(range(11), repeat=2))
    test_set = ((0,0),) + tuple((u,0) for u in range(1,5)) + tuple((0,u) for u in range(1,5))
    groups.append(('C11xC11_9_point_domain', elems, lambda x,y: ((x[0]+y[0])%11,(x[1]+y[1])%11), (0,0), test_set))
    source_total = candidate_total = 0
    for name, elems, mul, one, restricted in groups:
        powers, orders, cyclic = group_data(elems, mul, one)
        domain = tuple(range(len(elems))) if restricted is None else tuple(elems.index(x) for x in restricted)
        counts_source, counts_candidate = Counter(), Counter()
        limit = 6 if restricted is None else 9
        for n in range(1, min(limit, len(domain))+1):
            if any(2 <= o <= n+1 for o in orders):
                continue
            for a in combinations(domain, n):
                bounded = max(Counter(cyclic[x] for x in a).values()) <= 4
                assert has_assignment(a, powers), (name, n, a)
                if n <= 6:
                    counts_source[n] += 1
                if bounded:
                    counts_candidate[n] += 1
        source_total += sum(counts_source.values())
        candidate_total += sum(counts_candidate.values())
        print(f'GROUP {name} order={len(elems)} domain={len(domain)} '
              f'source_n<=6={dict(sorted(counts_source.items()))} '
              f'candidate={dict(sorted(counts_candidate.items()))}')
    print(f'GROUP_TOTAL source_fidelity_cases={source_total} candidate_cases={candidate_total}')


def main():
    faithful_checks()
    classification_checks()
    local_checks()
    group_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
