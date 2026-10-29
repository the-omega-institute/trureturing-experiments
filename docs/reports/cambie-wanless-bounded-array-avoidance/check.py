#!/usr/bin/env python3
"""Exact finite checks for Cambie--Wanless, arXiv:2610.20008v1, Conjecture 5.

No floating point, external dependencies, randomness, or claimed asymptotic proof.
Run: python3 check.py
"""
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial


def latin_squares(n):
    rows = tuple(permutations(range(n)))
    return tuple(L for L in product(rows, repeat=n)
                 if all(len({L[i][j] for i in range(n)}) == n for j in range(n)))


def principal_orbit(L):
    n = len(L)
    return {tuple(L[rho[i]][tau[j]] for i in range(n) for j in range(n))
            for rho in permutations(range(n)) for tau in permutations(range(n))}


def avoids(flat, A):
    return all(not (A[k] >> x) & 1 for k, x in enumerate(flat))


def bounded(A, n, t):
    return (all(x.bit_count() <= t for x in A)
            and all(sum((A[i*n+j] >> s) & 1 for j in range(n)) <= t
                    for i in range(n) for s in range(n))
            and all(sum((A[i*n+j] >> s) & 1 for i in range(n)) <= t
                    for j in range(n) for s in range(n)))


def conjugate_array(A, n):
    return tuple(sum(1 << j for j in range(n) if (A[i*n+j] >> s) & 1)
                 for i in range(n) for s in range(n))


def conjugate_square(L):
    n = len(L)
    return tuple(tuple(L[i].index(s) for s in range(n)) for i in range(n))


def partial_arrays(n):
    """All (1,1,1)-arrays, by extending one cell without row/symbol collisions."""
    A, row, col = [0]*(n*n), [0]*n, [0]*n
    def rec(k):
        if k == n*n:
            yield tuple(A)
            return
        i, j = divmod(k, n)
        A[k] = 0
        yield from rec(k+1)
        for s in range(n):
            bit = 1 << s
            if not ((row[i] | col[j]) & bit):
                A[k] = bit
                row[i] |= bit
                col[j] |= bit
                yield from rec(k+1)
                row[i] ^= bit
                col[j] ^= bit
        A[k] = 0
    yield from rec(0)


def source_fidelity():
    squares = {n: latin_squares(n) for n in range(1, 4)}
    assert [len(squares[n]) for n in range(1, 4)] == [1, 2, 12]
    # Section 4, the displayed partial Latin square, with symbols shifted by -1.
    P = (1, 2, 0, 2, 1, 4, 0, 4, 1)
    assert bounded(P, 3, 1)
    assert sum(avoids(tuple(x for row in L for x in row), P)
               for L in squares[3]) == 0
    print('FIDELITY: displayed order-3 partial Latin square: 0 of 12 Latin squares avoid it')
    pairs = 0
    for n in (2, 3):
        for L in squares[n]:
            orbit = principal_orbit(L)
            for Q in squares[n]:
                A = tuple(1 << x for row in Q for x in row)
                assert any(avoids(T, A) for T in orbit)
                pairs += 1
    print(f'FIDELITY: source Conjecture-4 small-order assertion checked on all {pairs} ordered pairs at n=2,3')
    return squares


def rook_checks():
    total = 0
    for n in range(4):
        perms = tuple(permutations(range(n)))
        for mask in range(1 << (n*n)):
            edges = tuple((i, j) for i in range(n) for j in range(n)
                          if (mask >> (i*n+j)) & 1)
            m = len(edges)
            D = max([0] + [sum(i == u for i, j in edges) for u in range(n)]
                    + [sum(j == v for i, j in edges) for v in range(n)])
            q = [sum(len({i for i, j in E}) == k
                     and len({j for i, j in E}) == k
                     for E in combinations(edges, k)) for k in range(n+1)]
            exact = Fraction(sum(all((i, pi[i]) not in edges for i in range(n))
                                 for pi in perms), factorial(n))
            prefix = Fraction(0)
            for k in range(n+1):
                prefix += (-1)**k * Fraction(q[k]*factorial(n-k), factorial(n))
                assert prefix <= exact if k % 2 else exact <= prefix
                if k >= 2 and m:
                    # Count ordered edge tuples with a conflicting pair.
                    defect = m**k - factorial(k)*q[k]
                    assert 0 <= defect <= comb(k, 2)*(2*D-1)*m**(k-1)
            assert prefix == exact
            total += 1
    print(f'ROOK: exact inclusion-exclusion, every Bonferroni truncation and collision bound: all {total} graphs at n=0..3')


def overlap_checks():
    n, r, checked = 2, 3, 0
    graphs = [frozenset(k for k in range(n*n) if (mask >> k) & 1)
              for mask in range(1 << (n*n))]
    def degree(F):
        return max([0] + [sum(k//n == i for k in F) for i in range(n)]
                   + [sum(k % n == j for k in F) for j in range(n)])
    for t in (1, 2):
        eligible = [F for F in graphs if degree(F) <= t]
        for family in product(eligible, repeat=r):
            d = [sum(k in F for F in family) for k in range(n*n)]
            if max(d) > t:
                continue
            lhs = sum(len(F & G) for F, G in combinations(family, 2))
            assert lhs == sum(comb(x, 2) for x in d)
            assert 2*lhs <= (t-1)*sum(map(len, family)) <= r*t*(t-1)*n
            checked += 1
    print(f'OVERLAP: double-counting and O(r) overlap bound: all {checked} admissible 3-graph families at n=2, t=1,2')


def hall_checks():
    tested, qualifying = 0, 0
    for n in range(1, 5):
        full = (1 << n)-1
        perms = tuple(permutations(range(n)))
        for rows in product(range(1 << n), repeat=n):
            delta = max([x.bit_count() for x in rows]
                        + [sum((row >> j) & 1 for row in rows) for j in range(n)])
            common = {}
            for r in range(1, n+1):
                sizes = []
                for X in combinations(range(n), r):
                    intersection = full
                    for x in X:
                        intersection &= rows[x]
                    sizes.append(intersection.bit_count())
                common[r] = max(sizes)
            # General sufficient criterion used in the proposed Hall argument.
            hypotheses = any(delta <= n-a and common[r] < a
                             for a in range(1, n+1) for r in range(1, a+1))
            if hypotheses:
                assert any(all(not (rows[i] >> pi[i]) & 1 for i in range(n))
                           for pi in perms)
                qualifying += 1
            tested += 1
    print(f'HALL: all {tested} bipartite graphs at n=1..4; {qualifying} satisfy the degree/common-neighborhood criterion and have complement matching')


def structural_checks(squares):
    counts = []
    for n in range(1, 4):
        # All arrays at n<=2; all (1,1,1)-arrays at n=3.
        arrays = tuple(product(range(1 << n), repeat=n*n)) if n <= 2 else tuple(partial_arrays(n))
        perms = tuple(permutations(range(n)))
        for A in arrays:
            Astar = conjugate_array(A, n)
            assert conjugate_array(Astar, n) == A
            for t in range(n+1):
                assert bounded(A, n, t) == bounded(Astar, n, t)
        for L in squares[n]:
            Lstar = conjugate_square(L)
            assert Lstar in squares[n]
            assert conjugate_square(Lstar) == L
            orbit = principal_orbit(L)
            # Invert the row-plus-symbol isotopes of Lstar, cell by cell.
            conjugated_orbit = set()
            for rho in perms:
                for sigma in perms:
                    T = tuple(tuple(sigma[Lstar[rho[i]][s]] for s in range(n))
                              for i in range(n))
                    conjugated_orbit.add(tuple(T[i].index(j) for i in range(n)
                                               for j in range(n)))
            assert conjugated_orbit == orbit
            for A in arrays:
                Astar = conjugate_array(A, n)
                # F(s,x) is independent of pi; construct it once per square/array.
                Fs = {(s, x): {(i, u) for i in range(n) for u in range(n)
                               if any((A[i*n+j] >> s) & 1 and L[u][j] == x
                                      for j in range(n))}
                      for s in range(n) for x in range(n)}
                for t in range(n+1):
                    if not bounded(A, n, t):
                        continue
                    for s in range(n):
                        for x in range(n):
                            F = Fs[s, x]
                            assert all(sum(i == v for i, u in F) <= t for v in range(n))
                            assert all(sum(u == v for i, u in F) <= t for v in range(n))
                        assert all(sum((i, u) in Fs[s, x] for x in range(n)) <= t
                                   for i in range(n) for u in range(n))
                for pi in perms:
                    # Direct cell-incidence graph, independent of the graph-hit test.
                    B = [0]*n
                    for i in range(n):
                        for j in range(n):
                            B[L[pi[i]][j]] |= A[i*n+j]
                    for s in range(n):
                        for x in range(n):
                            hit = any((i, pi[i]) in Fs[s, x] for i in range(n))
                            assert hit == bool((B[x] >> s) & 1)
                    # A complement matching is exactly an allowed symbol permutation.
                    for sigma in perms:
                        matching = all(not (B[x] >> sigma[x]) & 1 for x in range(n))
                        direct = all(not (A[i*n+j] >> sigma[L[pi[i]][j]]) & 1
                                     for i in range(n) for j in range(n))
                        assert matching == direct
        counts.append((n, len(squares[n]), len(arrays)))
    print(f'STRUCTURE: (order, all Latin squares, arrays)={counts}; conflict/matching identities, graph degrees, edge multiplicity and coordinate conjugation all exact')


def main():
    squares = source_fidelity()
    rook_checks()
    overlap_checks()
    hall_checks()
    structural_checks(squares)
    print('SCOPE: these finite checks test proof identities; they do not certify a uniform asymptotic threshold N(t)')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
