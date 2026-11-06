#!/usr/bin/env python3
"""Exact checks for Dement v5, Remark 7.8; no numerical approximations.

The reflection enumeration tests all fixed vertices and exchanged pairs on
all three axes using the three squared distances. It does not use the
proposed axis reduction or automaton to select triangles. For n <= 3 a
second scan of all unordered triples verifies completeness independently.
"""
from itertools import product, combinations
from fractions import Fraction

LAMBDA = ((0, 0), (1, 0), (0, 1), (-1, -1))
GAMMA = (0, 2, 3, 1)
TRANSPOSITIONS = ((0, 1, 3, 2), (0, 3, 2, 1), (0, 2, 1, 3))
TABLE3 = (
    (1, 1, 0, 1, 1, 0),
    (18, 9, 9, 6, 6, 0),
    (106, 43, 63, 25, 34, 9),
    (498, 165, 333, 87, 138, 51),
    (2146, 571, 1575, 274, 490, 216),
    (8898, 1869, 7029, 813, 1626, 813),
)

def centroid(w):
    """The source sum, evaluated directly (not the carry recursion)."""
    n = len(w)
    x = y = 0
    for r, digit in enumerate(w):
        sign = (-1) ** w[:r].count(0)
        weight = sign * 2 ** (n - r - 1)
        lx, ly = LAMBDA[digit]
        x += weight * lx
        y += weight * ly
    return x, y


def qdist(a, b):
    x, y = a[0] - b[0], a[1] - b[1]
    return x * x - x * y + y * y


def equilateral(T, coordinates):
    a, b, c = [coordinates[w] for w in T]
    k = qdist(a, b)
    return k > 0 and k == qdist(a, c) == qdist(b, c)


def local(T):
    # Literal synchronized gamma_S orbit, rather than a product parity test.
    a, b, c = T
    support = [r for r in range(len(a)) if a[r] != b[r]]
    if not support or any(a[r] == 0 for r in support):
        return False
    g = tuple(GAMMA[d] if r in support else d for r, d in enumerate(a))
    gg = tuple(GAMMA[d] if r in support else d for r, d in enumerate(g))
    return {a, b, c} == {a, g, gg}


def centered(T, coordinates):
    a, b, c = T
    unsigned_product = tuple(x ^ y ^ z for x, y, z in zip(a, b, c))
    p = coordinates[unsigned_product]
    return all(3 * p[k] == sum(coordinates[w][k] for w in T) for k in (0, 1))


def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def reflection_triangles(words, coordinates):
    found = set()
    per_axis = []
    for tau in TRANSPOSITIONS:
        image = {w: tuple(tau[d] for d in w) for w in words}
        fixed = [a for a in words if image[a] == a]
        pairs = [(b, image[b]) for b in words if b < image[b]]
        axis = set()
        for a in fixed:
            for b, c in pairs:
                T = tuple(sorted((a, b, c)))
                if equilateral(T, coordinates):
                    axis.add(T)
        found.update(axis)
        per_axis.append(axis)
    return found, per_axis


def brute_checks():
    print('Exact source checks: e,i,j,k coded 0,1,2,3; unsigned product is XOR.')
    print('n  sym local nonlocal local_product product nonlocal_product per_axis_product')
    for n, row in enumerate(TABLE3, 1):
        words = list(product(range(4), repeat=n))
        coordinates = {w: centroid(w) for w in words}
        assert len(set(coordinates.values())) == 4 ** n
        for w in words:
            sign = (-1) ** w.count(0)
            assert (sum(coordinates[w]) - (2 ** n - (1 if sign == 1 else 2))) % 3 == 0
        triangles, per_axis = reflection_triangles(words, coordinates)
        if n <= 3:
            all_equilateral = {T for T in combinations(words, 3) if equilateral(T, coordinates)}
            assert len(all_equilateral) == 4 ** n * (4 ** n - 1) // 12
            def symmetric(T):
                return any({tuple(tau[d] for d in w) for w in T} == set(T)
                           for tau in TRANSPOSITIONS)
            assert triangles == {T for T in all_equilateral if symmetric(T)}
            assert sum(map(local, all_equilateral)) == (7 ** n - 4 ** n) // 3
        loc = {T for T in triangles if local(T)}
        pc = {T for T in triangles if centered(T, coordinates)}
        actual = (len(triangles), len(loc), len(triangles - loc),
                  len(loc & pc), len(pc), len(pc - loc))
        assert actual == row, (n, actual, row)
        axis_pc = [sum(centered(T, coordinates) for T in axis) for axis in per_axis]
        assert len(set(axis_pc)) == 1
        assert len(triangles) == 9 * 4 ** (n - 1) - 5 * 2 ** n + 2
        assert len(loc & pc) == 3 * fib(2 * n + 2) - 5 * 2 ** n + 2
        if n == 1:
            assert (len(pc), len(pc - loc)) == (1, 0)
        else:
            assert axis_pc[0] == 8 * 3 ** (n - 2) - 2 ** n
            assert len(pc) == 8 * 3 ** (n - 1) - 5 * 2 ** n + 2
            assert len(pc - loc) == 8 * 3 ** (n - 1) - 3 * fib(2 * n + 2)
        print(n, *actual, axis_pc[0])
    print('ALL_SIX_TABLE3_COLUMNS_MATCH_n_1_TO_6')
    print('ALL_UNORDERED_TRIPLES_CROSS_CHECK_n_1_TO_3')


STATES = ((-1,-1,-1,1,-1), (-1,-1,-1,1,1),
          (0,0,-1,-1,-1), (0,0,-1,-1,1),
          (0,0,1,1,-1), (0,0,1,1,1),
          (1,1,1,-1,-1), (1,1,1,-1,1))
EDGE_ROWS = (
    ((1,3,1),),
    ((0,3,5),(1,3,0)),
    ((0,0,5),(0,1,7),(1,0,0),(1,1,2),(1,2,3)),
    ((0,0,4),(0,1,6),(0,2,7),(1,2,2)),
    ((0,0,3),(0,1,1),(0,2,0),(1,2,5)),
    ((0,0,2),(0,1,0),(1,0,7),(1,1,5),(1,2,4)),
    ((0,3,2),(1,3,7)),
    ((1,3,6),),
)


def transition(state, alpha, beta):
    d1, d2, sa, sb, sp = state
    pi = alpha ^ int(beta in (2, 3))
    lx, ly = LAMBDA[beta]
    return (2*d1 + sa*alpha - sb*(lx+ly),
            2*d2 + sp*pi - sb*lx,
            sa*(-1)**(alpha == 0), sb*(-1)**(beta == 0), sp*(-1)**(pi == 0))


def certificate_checks():
    expected = {(i, a, b, j) for i, row in enumerate(EDGE_ROWS) for a, b, j in row}
    assert len(expected) == 24
    found = set()
    for j, (t1,t2,A,B,C) in enumerate(STATES):
        for alpha, beta in product(range(2), range(4)):
            pi = alpha ^ int(beta in (2,3))
            lx, ly = LAMBDA[beta]
            sa, sb, sp = A*(-1)**(alpha == 0), B*(-1)**(beta == 0), C*(-1)**(pi == 0)
            d1 = Fraction(t1 - sa*alpha + sb*(lx+ly), 2)
            d2 = Fraction(t2 - sp*pi + sb*lx, 2)
            if d1.denominator != 1 or d2.denominator != 1:
                continue
            d1, d2 = int(d1), int(d2)
            if (d1 - d2) % 2:
                continue
            pre = (d1,d2,sa,sb,sp)
            assert pre in STATES, (j,alpha,beta,pre)
            i = STATES.index(pre)
            assert transition(pre,alpha,beta) == STATES[j]
            found.add((i,alpha,beta,j))
    assert found == expected
    print('INVERSE_TRANSITIONS: 64 target/letter cases, exactly 24 integral parity-compatible edges')
    M = [[0]*8 for _ in range(8)]
    for i,a,b,j in found:
        M[i][j] += 1
    def mv(v):
        return [sum(m*x for m,x in zip(row,v)) for row in M]
    w = [0,0,1,1,1,1,0,0]
    v1, v2 = mv(w), mv(mv(w))
    v3 = mv(v2)
    assert v1 == [0,1,3,2,2,3,1,0]
    assert v2 == [1,3,8,6,6,8,3,1]
    assert v3 == [3*x for x in v2]
    print('M*w =',v1)
    print('M^2*w =',v2)
    print('M^3*w =',v3,'= 3*M^2*w')
    # An independent unrestricted integer-state DP checks that no accepted
    # runs were lost. It applies no carry bound or pruning.
    distribution = {STATES[5]: 1}
    matrix_vector = w
    for n in range(1,11):
        next_distribution = {}
        for state, multiplicity in distribution.items():
            for alpha,beta in product(range(2),range(4)):
                t = transition(state,alpha,beta)
                assert (t[0]-t[1]) % 2 == 0
                next_distribution[t] = next_distribution.get(t,0) + multiplicity
        distribution = next_distribution
        accepted = sum(c for s,c in distribution.items() if s[0] == s[1] == 0)
        assert all(s[2] == s[3] for s in distribution if s[0] == s[1] == 0)
        matrix_vector = mv(matrix_vector)
        assert accepted == matrix_vector[5]
        assert accepted == (3 if n == 1 else 8 * 3 ** (n-2))
    print('UNPRUNED_INTEGER_STATE_DP_AND_MATRIX_MATCH_n_1_TO_10')


if __name__ == '__main__':
    brute_checks()
    certificate_checks()
    print('ALL_CHECKS_OK')
