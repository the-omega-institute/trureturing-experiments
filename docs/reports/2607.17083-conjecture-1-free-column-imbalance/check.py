#!/usr/bin/env python3
"""Exact checks for Serra et al., arXiv:2607.17083v1, Conjecture 1.

Python 3 and z3-solver 4.15.4. All mathematical data and solver constraints
are integers/Booleans. No floating-point optimization or tolerance is used.
"""
from itertools import combinations, product
from math import comb
import z3


def interactions(rows):
    k = len(rows[0])
    return {(i, j, r[i], r[j]) for r in rows for i, j in combinations(range(k), 2)}


def constraints(b):
    return {(i, b + j, 1, 1) for i in range(b) for j in range(b) if i != j}


def valid_rows(b):
    forbidden = constraints(b)
    return [r for r in product(range(2), repeat=2 * b + 1)
            if all((r[i], r[j]) != (a, c) for i, j, a, c in forbidden)]


def required_pairs(b):
    return {(i, j, a, c) for i, j in combinations(range(2 * b + 1), 2)
            for a, c in product(range(2), repeat=2)} - constraints(b)


def displayed(b):
    a = (1,) * b + (0,) * b + (1,)
    c = (0,) * b + (1,) * b + (1,)
    return [a, c] + [tuple(int(j in (i, b + i)) for j in range(2 * b)) + (0,)
                     for i in range(b)]


def covers(b, rows):
    return (all(all((r[i], r[j]) != (a, c) for i, j, a, c in constraints(b))
                for r in rows) and interactions(rows) == required_pairs(b))


def reducedness(b, rows):
    k = 2 * b + 1
    assert len(rows) == 2 * (2 ** (b + 1) + b - 1)
    for i in range(k):
        assert {r[i] for r in rows} == {0, 1}
    # Enumerate complete valid assignments, not only the proposed cover.
    assert interactions(rows) == required_pairs(b)
    assert all(sum((i, j, a, c) in constraints(b) for a, c in product(range(2), repeat=2)) <= 1
               for i, j in combinations(range(k), 2))


def integer_cover_solver(b):
    rows = valid_rows(b)
    reducedness(b, rows)
    counts = [z3.Int(f'row_{n}') for n in range(len(rows))]
    s = z3.SolverFor('QF_LIA')
    s.add(*(v >= 0 for v in counts))
    for i, j, a, c in sorted(required_pairs(b)):
        s.add(z3.Sum([v for v, r in zip(counts, rows) if (r[i], r[j]) == (a, c)]) >= 1)
    total = z3.Sum(counts)
    ones = z3.Sum([v for v, r in zip(counts, rows) if r[-1] == 1])
    return rows, counts, s, total, ones


def expect(s, condition, expected):
    s.push()
    s.add(condition)
    result = s.check()
    assert result == expected, (str(condition), result, expected)
    witness = s.model() if result == z3.sat else None
    s.pop()
    return witness


def fidelity():
    # Section 2, displayed matrix M^(1), copied by rows from the TeX.
    matrix = [
        (1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
        (0, 0, 0, 1, 0, 0, 1, 0, 1, 1),
        (0, 0, 1, 0, 0, 1, 0, 1, 0, 1),
        (0, 1, 0, 0, 1, 0, 0, 1, 1, 0),
        (1, 0, 0, 0, 1, 1, 1, 0, 0, 0),
        (1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
    ]
    assert len(interactions(matrix)) == 4 * comb(10, 2) == 180
    assert all(sum(r[j] for r in matrix) == 3 for j in range(10))
    assert all(comb(n - 1, (n + 1) // 2) < 10 for n in range(1, 6))
    assert comb(5, 3) == 10
    print('FIDELITY source M^(1): 6 rows, 10 columns, 180 interactions, CAN(2,10,2)=6')


def family_checks():
    for b in range(3, 10):
        candidate = displayed(b)
        assert covers(b, candidate)
        assert sorted([sum(r[-1] == 0 for r in candidate), sum(r[-1] == 1 for r in candidate)]) == [2, b]
        print(f'FAMILY b={b}: {2*b+1} parameters, {b*(b-1)} constraints, '
              f'{b+2} valid covering rows, free-column frequencies=(2,{b})')


def exact_search(b):
    rows, counts, s, total, ones = integer_cover_solver(b)
    reducedness(b, rows)
    optimum = b + 2
    expect(s, total <= optimum - 1, z3.unsat)
    expect(s, total == optimum, z3.sat)
    excluded = []
    for w in range(optimum + 1):
        result = z3.sat if w in (2, b) else z3.unsat
        expect(s, z3.And(total == optimum, ones == w), result)
        if result == z3.unsat:
            excluded.append(w)
    # Enumerate all optimal integer-multiplicity vectors, including any
    # putative optimum containing duplicate rows.
    s.push()
    s.add(total == optimum)
    solutions = []
    while s.check() == z3.sat:
        m = s.model()
        values = [m.eval(v).as_long() for v in counts]
        assert all(v in (0, 1) for v in values)
        chosen = frozenset(r for v, r in zip(values, rows) if v)
        assert covers(b, list(chosen))
        solutions.append(chosen)
        s.add(z3.Or([v != n for v, n in zip(counts, values)]))
    s.pop()
    standard = frozenset(displayed(b))
    complemented = frozenset(r[:-1] + (1-r[-1],) for r in standard)
    assert set(solutions) == {standard, complemented}
    print(f'EXACT b={b}: valid complete tests={len(rows)}, required interactions={len(required_pairs(b))}; '
          f'rows<={optimum-1}: UNSAT; optimum={optimum}; '
          f'all optimal row-multisets={len(solutions)}; possible z-one-counts={sorted({2,b})}; '
          f'other counts {excluded}: UNSAT')
    if b == 4:
        for n in range(1, 7):
            expect(s, z3.And(total == n, 2*ones >= n-1, 2*ones <= n+1), z3.unsat)
        m = expect(s, z3.And(total == 7, 2*ones >= 6, 2*ones <= 8), z3.sat)
        assert m is not None
        augmented = displayed(4) + [(0,)*8 + (1,)]
        assert covers(4, augmented)
        assert sum(r[-1] for r in augmented) == 3
        print('BALANCE b=4: every balanced cover with <=6 rows UNSAT; 7 rows SAT; explicit z counts=(4,3)')
        print('MATRIX b=4 (x_1,...,x_4,y_1,...,y_4,z):')
        for r in displayed(4):
            print(' '.join(map(str, r)))


def small_boundaries():
    for b, expected in [(1, 4), (2, 5)]:
        rows, counts, s, total, ones = integer_cover_solver(b)
        expect(s, total < expected, z3.unsat)
        expect(s, z3.And(total == expected, 2*ones >= expected-1,
                         2*ones <= expected+1), z3.sat)
        assert not covers(b, displayed(b))
        print(f'BOUNDARY b={b}: optimum={expected}; balanced optimum SAT; displayed b+2 rows fail coverage')
    print('BOUNDARY b=3: optimum=5; both optimal z columns have counts=(2,3), hence balanced')


def main():
    print('Arithmetic: integer/Boolean only; Z3 ' + z3.get_version_string() + ' QF_LIA')
    fidelity()
    family_checks()
    for b in (3, 4, 5):
        exact_search(b)
    small_boundaries()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
