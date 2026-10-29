#!/usr/bin/env python3
"""Exact checks for the direct-sum clause of KMV Question 12.

No floating point. Finite checks support, but do not prove, the all-field family.
Source fidelity: KMV arXiv:2610.22248v1, Section 1.1, first example:
if e_i^2=e_i, the function algebra K^n has 2^n idempotents.
"""
from fractions import Fraction
from itertools import combinations, product
import sympy as sp


def rank(rows, modulus=None):
    if not rows:
        return 0
    a = [[Fraction(x) if modulus is None else int(x) % modulus for x in r]
         for r in rows]
    height, width, pivot = len(a), len(a[0]), 0
    for j in range(width):
        k = next((i for i in range(pivot, height) if a[i][j] != 0), None)
        if k is None:
            continue
        a[pivot], a[k] = a[k], a[pivot]
        inv = 1 / a[pivot][j] if modulus is None else pow(a[pivot][j], -1, modulus)
        a[pivot] = [x * inv if modulus is None else x * inv % modulus for x in a[pivot]]
        for i in range(height):
            if i == pivot:
                continue
            c = a[i][j]
            a[i] = [x - c * y if modulus is None else (x - c * y) % modulus
                    for x, y in zip(a[i], a[pivot])]
        pivot += 1
        if pivot == height:
            break
    return pivot


def unit(d, i):
    return tuple(int(j == i) for j in range(d))


def add(x, y, modulus=None):
    z = tuple(a + b for a, b in zip(x, y))
    return z if modulus is None else tuple(a % modulus for a in z)


def natural_product(squares, x, y, modulus=None):
    z = tuple(sum(a * b * q[j] for a, b, q in zip(x, y, squares))
              for j in range(len(x)))
    return z if modulus is None else tuple(a % modulus for a in z)


def factor(d, negative_fourth=False):
    w = add(unit(d, 0), unit(d, 1))
    signs = [1, -1] + [1] * (d - 2)
    if negative_fourth:
        signs[3] = -1
    return [tuple(s * a for a in w) for s in signs], signs


def annihilator_rank(squares, modulus=None):
    d = len(squares)
    # Row block j is multiplication by the j-th natural basis vector.
    rows = [[squares[j][k] if i == j else 0 for i in range(d)]
            for j in range(d) for k in range(d)]
    return rank(rows, modulus)


def fidelity():
    cases = 0
    for p in (3, 5):
        for n in range(1, 6):
            squares = [unit(n, i) for i in range(n)]
            idem = [v for v in product(range(p), repeat=n)
                    if natural_product(squares, v, v, p) == v]
            assert len(idem) == 2 ** n
            assert set(idem) == set(product((0, 1), repeat=n))
            cases += 1
    print(f"FIDELITY: function-algebra idempotents = 2^n; F_3,F_5; n=1..5; {cases} cases")


def family():
    cases = 0
    table = {(0, 0): (1, 1), (0, 1): (0, 1),
             (1, 0): (0, 1), (1, 1): (-1, 1)}
    for m in range(4, 10):
        for n in range(3, 9):
            a, _ = factor(m, True)
            b, _ = factor(n)
            ambient = [q + (0,) * n for q in a] + [(0,) * m + q for q in b]
            x = add(unit(m + n, 2), unit(m + n, m + 2))
            y = add(unit(m + n, 3), unit(m + n, m + 2))
            z = add(unit(m + n, 0), unit(m + n, 1))
            t = add(unit(m + n, m), unit(m + n, m + 1))
            basis = [x, y, z, t]
            for p in (None, 3, 5, 7):
                assert rank(a, p) == rank(b, p) == 1
                assert annihilator_rank(a, p) == m
                assert annihilator_rank(b, p) == n
                assert rank(basis, p) == 4
                outputs = []
                for i, j in product(range(4), repeat=2):
                    c, d = table.get((i, j), (0, 0))
                    expected = tuple(c * zi + d * ti for zi, ti in zip(z, t))
                    if p is not None:
                        expected = tuple(v % p for v in expected)
                    actual = natural_product(ambient, basis[i], basis[j], p)
                    assert actual == expected, (m, n, p, i, j)
                    outputs.append(actual)
                assert rank(outputs, p) == 2
                cases += 1
    print(f"FAMILY: m=4..9,n=3..8 over Q,F_3,F_5,F_7; {cases} cases; all 16 ordered S-products each")
    print("FAMILY: dim A^2=dim B^2=1; Ann(A)=Ann(B)=0; dim S=4; dim S^2=2")


def subspaces(d, p):
    """Every subspace once, via its unique reduced-row-echelon basis."""
    for r in range(d + 1):
        for pivots in combinations(range(d), r):
            slots = [(i, j) for i, pivot in enumerate(pivots)
                     for j in range(pivot + 1, d) if j not in pivots]
            for values in product(range(p), repeat=len(slots)):
                rows = [list(unit(d, pivot)) for pivot in pivots]
                for (i, j), value in zip(slots, values):
                    rows[i][j] = value
                yield rows, pivots


def contained(v, rows, pivots, p):
    return all((v[j] - sum(v[pivot] * row[j]
                          for pivot, row in zip(pivots, rows))) % p == 0
               for j in range(len(v)))


def orthogonal_basis(rows, signs, p):
    def dot(v, w):
        return sum(s * x * y for s, x, y in zip(signs, v, w)) % p
    todo, answer = [r[:] for r in rows], []
    while todo:
        i = next((i for i, v in enumerate(todo) if dot(v, v)), None)
        if i is None:
            pair = next(((i, j) for i, j in combinations(range(len(todo)), 2)
                         if dot(todo[i], todo[j])), None)
            if pair is None:
                answer.extend(todo)
                break
            i, j = pair
            todo[i] = [(x + y) % p for x, y in zip(todo[i], todo[j])]
            assert dot(todo[i], todo[i]) != 0  # 2 is invertible.
        v = todo.pop(i)
        inv = pow(dot(v, v), -1, p)
        todo = [[(a - dot(v, w) * inv * b) % p for a, b in zip(w, v)]
                for w in todo]
        answer.append(v)
    assert rank(answer, p) == len(rows)
    assert all(dot(v, w) == 0 for i, v in enumerate(answer) for w in answer[i + 1:])
    return answer


def exhaustive_heredity():
    for p in (3, 5):
        for name, d, fourth in (('A', 4, True), ('B', 3, False)):
            squares, signs = factor(d, fourth)
            seen, closed = 0, 0
            for rows, pivots in subspaces(d, p):
                seen += 1
                if not all(contained(natural_product(squares, v, w, p), rows, pivots, p)
                           for v in rows for w in rows):
                    continue
                closed += 1
                natural = orthogonal_basis(rows, signs, p)
                assert all(contained(v, rows, pivots, p) for v in natural)
                assert all(natural_product(squares, v, w, p) == (0,) * d
                           for i, v in enumerate(natural) for w in natural[i + 1:])
            # Gaussian-binomial counts, evaluated independently for these dimensions.
            expected = {3: {3: 28, 4: 212}, 5: {3: 64, 4: 1120}}[p][d]
            assert seen == expected
            print(f"HEREDITY: {name}_{d} over F_{p}; {seen} subspaces, {closed} closed; all have verified natural bases")


def obstruction():
    a, b, c, d, r, s, h, k = sp.symbols('a b c d r s h k')
    aa, _ = factor(4, True)
    bb, _ = factor(3)
    squares = [q + (0,) * 3 for q in aa] + [(0,) * 4 + q for q in bb]
    left, right = (r, r, a, b, s, s, a + b), (h, h, c, d, k, k, c + d)
    out = natural_product(squares, left, right)
    f, g = a * c - b * d, (a + b) * (c + d)
    assert all(sp.expand(x - y) == 0 for x, y in zip(out, (f, f, 0, 0, g, g, 0)))
    q1, q2 = sp.diag(1, -1), sp.ones(2)
    operator = q1.inv() * q2
    assert operator != sp.zeros(2) and operator ** 2 == sp.zeros(2)
    print("SYMBOLIC: lifts multiply to (ac-bd)z+(a+b)(c+d)t; all annihilator coordinates cancel")
    print("SYMBOLIC: Q1^(-1)Q2=[[1,1],[-1,-1]] is nonzero and square-zero")
    for p in (3, 5, 7, 11, 13):
        zero_pairs = independent = 0
        for a, b, c, d in product(range(p), repeat=4):
            if (a * c - b * d) % p == 0 and ((a + b) * (c + d)) % p == 0:
                zero_pairs += 1
                if (a * d - b * c) % p != 0:
                    independent += 1
        assert independent == 0
        print(f"OBSTRUCTION: F_{p}; {p**4} ordered quotient pairs; {zero_pairs} zero-product pairs; 0 independent")


def main():
    print('Arithmetic: Python integers/Fraction and SymPy exact expressions; no floating point')
    fidelity()
    family()
    exhaustive_heredity()
    obstruction()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
