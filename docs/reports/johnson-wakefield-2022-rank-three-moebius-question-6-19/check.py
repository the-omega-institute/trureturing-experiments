"""Johnson-Wakefield, Ars Math. Contemp. 24 (2024) #P2.03 (arXiv:2204.07158v1), Question 6.19.

The J-Moebius polynomial of a finite lattice of flats L is M(L, t) = sum over flags x <= y <= z of J(x, y, z) t^(3 rk L
- rk x - rk y - rk z) (Definition 6.2). The script computes J from its defining system (Definition 5.1) by the
triangular recursion, independently of Theorem 5.2, and checks J = mu(x, y) mu(y, z). It reproduces the two
factorizations of Remark 6.18 (U_{3,4} and the dual of the graphic matroid of K_{3,3}), and checks M(L, -1) = 2(b - v)
on rank-3 matroids (v points, b lines), with M(L, -1) = 0 exactly for the modular ones.
"""
import itertools
import sys

import sympy as sp

t = sp.symbols('t')


def flats_of(E, rank):
    E = list(E)
    out = []
    for r in range(len(E) + 1):
        for S in itertools.combinations(E, r):
            S = frozenset(S)
            if all(rank(S | {e}) > rank(S) for e in E if e not in S):
                out.append(S)
    return out


def mobius(F):
    order = sorted(F, key=len)
    mu = {}
    for x in F:
        mu[(x, x)] = 1
        for y in order:
            if x < y:
                mu[(x, y)] = -sum(mu[(x, z)] for z in order if x <= z < y)
    return mu


def J_by_definition(F):
    """Definition 5.1: sum over x <= a <= y <= b <= z of J(a, y, b) = delta_3(x, y, z), solved by recursion."""
    J = {}
    flags = [(x, y, z) for x in F for y in F for z in F if x <= y <= z]
    flags.sort(key=lambda f: len(f[2]) - len(f[0]))
    for (x, y, z) in flags:
        s = sum(J[(a, y, b)] for a in F for b in F if x <= a <= y <= b <= z and (a, b) != (x, z))
        J[(x, y, z)] = (1 if x == y == z else 0) - s
    return J


def Mpoly(E, rank):
    F = flats_of(E, rank)
    mu = mobius(F)
    J = J_by_definition(F)
    assert all(J[(x, y, z)] == mu[(x, y)] * mu[(y, z)] for (x, y, z) in J)
    R = rank(frozenset(E))
    P = sum(J[(x, y, z)] * t ** (3 * R - rank(x) - rank(y) - rank(z)) for (x, y, z) in J)
    return sp.expand(P), F


def uniform(r, n):
    return list(range(n)), (lambda S: min(len(S), r))


def rank3_from_lines(n, lines):
    L = [frozenset(l) for l in lines]

    def rank(S):
        S = frozenset(S)
        if len(S) <= 2:
            return len(S)
        return 2 if any(S <= l for l in L) else 3
    return list(range(n)), rank


def graphic_dual(nv, edges):
    def r(S):
        parent = list(range(nv))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        c = 0
        for i in S:
            u, w = edges[i]
            a, b = find(u), find(w)
            if a != b:
                parent[a] = b
                c += 1
        return c
    E = list(range(len(edges)))
    rE = r(E)
    return E, (lambda S: len(S) + r([e for e in E if e not in S]) - rE)


ok = True
E, r = uniform(3, 4)
P, _ = Mpoly(E, r)
target = sp.expand((t - 1) * (t**8 - 3*t**7 - t**6 + 12*t**5 - 2*t**4 - 12*t**3 + 3*t**2 + 5*t - 1))
print("Remark 6.18, U_{3,4}:", sp.factor(P), "matches source:", P == target)
ok &= P == target
K33 = [(i, 3 + j) for i in range(3) for j in range(3)]
E, r = graphic_dual(6, K33)
P, _ = Mpoly(E, r)
target = sp.expand((t**10 - 9*t**9 + 22*t**8 + 12*t**7 - 81*t**6 + 21*t**5 + 69*t**4 - 18*t**3 - 34*t**2 + 15*t - 1)
                   * (t + 1) * (t - 1))
print("Remark 6.18, M*(K_{3,3}):", "matches source:", P == target)
ok &= P == target

examples = {
    "U_{3,4}": uniform(3, 4), "U_{3,5}": uniform(3, 5), "U_{3,6}": uniform(3, 6),
    "Fano": rank3_from_lines(7, [(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (1, 4, 6), (2, 3, 6), (2, 4, 5)]),
    "non-Fano": rank3_from_lines(7, [(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (1, 4, 6), (2, 3, 6)]),
    "near-pencil(5)": rank3_from_lines(5, [(0, 1, 2, 3)]),
    "P_7 triangle": rank3_from_lines(6, [(0, 1, 2), (2, 3, 4), (4, 5, 0)]),
    "two disjoint lines": rank3_from_lines(6, [(0, 1, 2), (3, 4, 5)]),
    "AG(2,3)": rank3_from_lines(9, [(0, 1, 2), (3, 4, 5), (6, 7, 8), (0, 3, 6), (1, 4, 7), (2, 5, 8),
                                    (0, 4, 8), (1, 5, 6), (2, 3, 7), (0, 5, 7), (1, 3, 8), (2, 4, 6)]),
}
for name, (E, r) in examples.items():
    P, F = Mpoly(E, r)
    v = sum(1 for f in F if r(f) == 1)
    b = sum(1 for f in F if r(f) == 2)
    lines = [f for f in F if r(f) == 2]
    modular = all(len(a & c) > 0 for a, c in itertools.combinations(lines, 2))
    val = P.subs(t, -1)
    good = val == 2 * (b - v) and ((val == 0) == modular)
    ok &= good
    print(f"{name}: v={v} b={b} M(L,-1)={val} 2(b-v)={2 * (b - v)} modular={modular} ok={good}")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
