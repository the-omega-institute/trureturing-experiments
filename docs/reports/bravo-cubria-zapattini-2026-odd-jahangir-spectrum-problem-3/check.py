#!/usr/bin/env python3
"""Exact finite checks for arXiv:2610.28861v1, Problem 3. No numerical spectra."""
from functools import lru_cache
import sympy as sp

x, q_symbol = sp.symbols('x q')


def require(test, label):
    if not test:
        raise AssertionError(label)


def eq(left, right, label):
    require(sp.Poly(left - right, x, q_symbol).is_zero, label)


def cycle(n):
    require(n >= 3, 'cycle length')
    a = sp.zeros(n)
    for j in range(n):
        a[j, (j + 1) % n] = a[(j + 1) % n, j] = 1
    return a


def jahangir(m, q):
    n = m * q
    a = sp.zeros(n + 1)
    a[1:, 1:] = cycle(n)
    for j in range(0, n, m):
        a[0, j + 1] = a[j + 1, 0] = 1
    return a


def quotient(a, m, q):
    # Construct classes directly from the graph, not the proposed quotient formula.
    classes = [[0]] + [[] for _ in range(m // 2 + 1)]
    for j in range(m * q):
        residue = j % m
        distance = min(residue, m - residue)
        classes[distance + 1].append(j + 1)
    b = sp.zeros(len(classes))
    for i, cell in enumerate(classes):
        require(bool(cell), 'nonempty quotient cell')
        for j, target in enumerate(classes):
            counts = {sum(a[v, w] for w in target) for v in cell}
            require(len(counts) == 1, 'equitable neighbour counts')
            b[i, j] = counts.pop()
    # Check the defining equitable-partition intertwining equation independently.
    p = sp.zeros(a.rows, len(classes))
    for j, cell in enumerate(classes):
        for v in cell:
            p[v, j] = 1
    require(a * p == p * b, 'A P = P Q')
    return b


def predicted_odd_quotient(r, q):
    b = sp.zeros(r + 2)
    b[0, 1], b[1, 0] = q, 1
    b[1, 2], b[2, 1] = 2, 1
    for j in range(2, r + 1):
        b[j, j + 1] = b[j + 1, j] = 1
    b[r + 1, r + 1] = 1
    return b


@lru_cache(None)
def path_poly(j):
    if j == -1:
        return sp.Integer(0)
    if j == 0:
        return sp.Integer(1)
    return sp.expand(x * path_poly(j - 1) - path_poly(j - 2))


@lru_cache(None)
def cycle_poly(n):
    return cycle(n).charpoly(x).as_expr()


def main():
    print('Arithmetic: SymPy exact integer/symbolic polynomial arithmetic; no floats.')
    # Source fidelity: section 4.5 m=2 formula; n >= 4.
    fidelity_two = 0
    for q in range(2, 9):
        a = jahangir(2, q)
        b = quotient(a, 2, q)
        pi = b.charpoly(x).as_expr()
        eq(pi, x**3 - (q + 4) * x, 'source m=2 quotient')
        eq(a.charpoly(x).as_expr() * (x**2 - 4), cycle_poly(2*q) * pi,
           'source m=2 factorization')
        fidelity_two += 1
    print(f'FIDELITY m=2, q=2..8: {fidelity_two} source factorizations OK')
    # Source fidelity: Corollary 4.6 (after Proposition 4.5), even m >= 4.
    fidelity_even = 0
    for m in (4, 6, 8, 10):
        for q in range(1, 7):
            a = jahangir(m, q)
            pi = quotient(a, m, q).charpoly(x).as_expr()
            inherited = path_poly(m // 2 - 1)
            eq(a.charpoly(x).as_expr() * cycle_poly(m),
               cycle_poly(m*q) * inherited * pi, 'source even factorization')
            fidelity_even += 1
    print(f'FIDELITY even m=4,6,8,10, q=1..6: {fidelity_even} factorizations OK')
    # Symbolic quotient formula in q; this is finite in r, not a uniform proof.
    for r in range(1, 11):
        ar = path_poly(r) + path_poly(r - 1)
        br = path_poly(r) - path_poly(r - 1)
        pi = predicted_odd_quotient(r, q_symbol).charpoly(x).as_expr()
        eq(pi, x * (x - 2) * ar - q_symbol * br, 'odd quotient formula')
        anti = sp.zeros(r)
        for j in range(r - 1):
            anti[j, j + 1] = anti[j + 1, j] = 1
        anti[r - 1, r - 1] = -1
        eq(anti.charpoly(x).as_expr(), ar, 'reflection minus block')
    print('SYMBOLIC r=1..10: quotient formula and reflection-minus block OK')
    # Exhaust the stated grid, including every q=1 boundary case.
    cases = 0
    for m in range(3, 16, 2):
        for q in range(1, 9):
            if m * q > 80:
                continue
            r = (m - 1) // 2
            a = jahangir(m, q)
            b = quotient(a, m, q)
            require(b == predicted_odd_quotient(r, q), 'odd quotient entries')
            ar = path_poly(r) + path_poly(r - 1)
            br = path_poly(r) - path_poly(r - 1)
            pi = b.charpoly(x).as_expr()
            eq(pi, x * (x - 2) * ar - q * br, 'odd pi')
            division, remainder = sp.div(cycle_poly(m*q), cycle_poly(m), x)
            eq(remainder, 0, 'cycle quotient is a polynomial')
            eq(a.charpoly(x).as_expr(), division * ar * pi, 'odd full identity')
            require(sp.degree(division, x) + sp.degree(ar, x) + sp.degree(pi, x)
                    == m*q + 1, 'multiplicity dimension count')
            eq(sp.gcd(sp.Poly(ar, x), sp.Poly(pi, x)).as_expr(), 1,
               'second inherited family coprime to quotient')
            print(f'ODD m={m} q={q} n={m*q}: identity/equitable quotient/degrees/gcd OK')
            cases += 1
    print(f'ODD_GRID_OK: {cases} cases, m=3,5,...,15; q=1..8; mq<=80')
    for r, expected in [
        (1, x**3-x**2-(q_symbol+2)*x+q_symbol),
        (2, x**4-x**3-(q_symbol+3)*x**2+(q_symbol+2)*x+q_symbol),
    ]:
        eq(predicted_odd_quotient(r, q_symbol).charpoly(x).as_expr(), expected,
           'displayed exact example')
    print('EXAMPLES m=3 and m=5: displayed symbolic formulas OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
