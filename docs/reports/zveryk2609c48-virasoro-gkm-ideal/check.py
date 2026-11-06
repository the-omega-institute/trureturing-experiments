#!/usr/bin/env python3
"""Exact finite checks for the proposed Zveryk Conjecture 4.8 counterexample.

Only fractions and integers are used. This does not establish the presentation's
triangular decomposition, ideal disjointness, or the infinite family theorem.
"""
from fractions import Fraction as F
from functools import lru_cache
from collections import defaultdict


def add(*parts):
    out = defaultdict(F)
    for a, v in parts:
        for w, c in v.items():
            out[w] += a * c
    return {w: c for w, c in out.items() if c}


@lru_cache(None)
def normal(word):
    """PBW order for negative modes: decreasing positive subscripts.

    L_-a L_-b = L_-b L_-a + (b-a)L_-(a+b) when a<b.
    """
    for i in range(len(word) - 1):
        a, b = word[i:i + 2]
        if a < b:
            swapped = word[:i] + (b, a) + word[i + 2:]
            merged = word[:i] + (a + b,) + word[i + 2:]
            return add((1, normal(swapped)), (b - a, normal(merged)))
    return {word: F(1)}


HIGHEST = -F(3, 8)
CENTRAL = F(4)


@lru_cache(None)
def mode(m, word):
    """Virasoro action on a PBW word in M(c=4,h=-3/8)."""
    if m < 0:
        return normal((-m,) + word)
    if m == 0:
        return {word: HIGHEST + sum(word)}
    if not word:
        return {}
    a, tail = word[0], word[1:]
    left = {}
    for w, c in mode(m, tail).items():
        left = add((1, left), (c, normal((a,) + w)))
    central = {tail: (F(m**3 - m, 12) * CENTRAL)} if m == a else {}
    return add((1, left), (m + a, mode(m - a, tail)), (1, central))


def act(m, v):
    return add(*[(c, mode(m, w)) for w, c in v.items()])


def partitions(n, cap=None):
    if n == 0:
        return [()]
    cap = min(n, n if cap is None else cap)
    return [(a,) + p for a in range(cap, 0, -1)
            for p in partitions(n - a, a)]


def rank(rows):
    a = [[F(x) for x in row] for row in rows]
    k = 0
    for j in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(k, len(a)) if a[i][j]), None)
        if pivot is None:
            continue
        a[k], a[pivot] = a[pivot], a[k]
        v = a[k][j]
        a[k] = [x / v for x in a[k]]
        for i in range(len(a)):
            if i != k:
                v = a[i][j]
                a[i] = [x - v * y for x, y in zip(a[i], a[k])]
        k += 1
    return k


def cartan_action(word, delta):
    """Iterated adjoint action on H; 0 denotes H, r>0 denotes L_-r."""
    v = {0: F(1)}
    for r in reversed(word):
        out = defaultdict(F)
        for k, c in v.items():
            if k == 0:
                out[r] += c * r * delta
            else:
                out[r + k] += c * (k - r)
        v = {k: c for k, c in out.items() if c}
    return v


def main():
    # Fidelity: source equation (3.1), r=s=2. The square-root term
    # has coefficient r^2-s^2=0 and is exactly absent, even at c=4.
    h22 = F((13 - 4) * (2**2 + 2**2) - 24*2*2 - 2 + 2*4, 48)
    assert h22 == HIGHEST
    # Fidelity: the source's presented bracket has 6C, hence c=12C.
    assert F(4) * HIGHEST + F(6) * CENTRAL / 12 == F(1, 2)
    assert mode(2, (2,)) == {(): F(1, 2)}
    print('FIDELITY: source h_(2,2)(4)=-3/8; [L_2,L_-2]=4B+6C; c=12C')

    for q in range(1, 101):
        a = ((-q-F(1,3), -F(q)), (-F(q), -q+F(1,3)))
        c = (F(1), -F(1))
        b = (-9*q-F(9,8), 9*q+F(33,8))
        gamma = (-F(1,3), -F(1,3))
        pair = lambda root, x: sum(root[i]*x[i] for i in range(2))
        roots = tuple(tuple(a[i][j] for i in range(2)) for j in range(2))
        assert a[0][0]*a[1][1]-a[0][1]*a[1][0] == -F(1,9)
        assert c[0]*b[1]-c[1]*b[0] == 3
        assert pair(gamma,c) == 0 and pair(gamma,b) == -1
        assert tuple(roots[0][i]-roots[1][i] for i in range(2)) == gamma
        assert [pair(r,c) for r in roots] == [-F(1,3), -F(1,3)]
        assert [pair(r,b) for r in roots] == [F(3,8), F(11,8)]
        assert -pair(roots[0],c)*12 == CENTRAL
        assert all(x < 0 for row in a for x in row)
        assert all(-r[0] > 0 for r in roots+(gamma,))
        for i in range(2):
            local_weights = [-2*a[i][1-i]/a[i][i]]
            local_weights += [F(2*m,3)/a[i][i] for m in range(1,9)]
            for lam in local_weights:
                assert lam < 0
                assert all(r*(lam-r+1) != 0 for r in range(1,10))
    print('FAMILY: q=1..100 exact Cartan, symmetry, (4.3), central/degree independence, and local-weight checks passed')

    p = {(1,1,1,1): F(1), (2,1,1): -F(3), (2,2): -F(7,4),
         (3,1): F(3), (4,): F(3,2)}
    assert all(sum(w) == 4 for w in p)
    assert p[(1,1,1,1)] == 1 and len(p) == 5
    for m in range(1,9):
        assert not act(m,p), (m,act(m,p))
    # Check the Virasoro action itself on the finite PBW window.
    count = 0
    for level in range(5):
        for w in partitions(level):
            for m in range(-2,3):
                for n in range(-2,3):
                    lhs = add((1,act(m,mode(n,w))),(-1,act(n,mode(m,w))))
                    rhs = add((m-n,mode(m+n,w)),
                              (F(m**3-m,12)*CENTRAL if m == -n else 0,{w:F(1)}))
                    assert lhs == rhs, (m,n,w,lhs,rhs)
                    count += 1
    cols = partitions(4)
    rows = [[mode(m,w).get(t,F(0)) for w in cols]
            for m in (1,2) for t in partitions(4-m)]
    assert rank(rows) == 4 and len(cols) == 5
    print('SINGULAR: q=1, five PBW coefficients; L_1..L_8 annihilate P f_1; joint L_1/L_2 kernel at level 4 has dimension 1')
    print(f'ACTION: {count} exact commutator checks on PBW levels 0..4, modes -2..2 passed')

    delta = -F(1,3)
    terms = [(coeff,cartan_action(w,delta)) for w,coeff in p.items()]
    assert not add(*terms)
    assert cartan_action((3,1),delta) == {4:-2*delta}
    assert cartan_action((4,),delta) == {4:4*delta}
    assert all(not cartan_action(w,delta) for w in [(1,1,1,1),(2,1,1),(2,2)])
    assert -6*delta+6*delta == 0
    print('CARTAN: P H_1=0, with final terms -6 delta L_-4 + 6 delta L_-4; e_2 commutes with every factor')
    # A negative-Witt detector needs no PBW independence theorem.
    # On Laurent monomials, L_-r(t^k)=-k*t^(k-r).
    # In the semidirect representation f_1 sends 1 to t and f_2 to zero.
    detector = F(0)
    for word, coeff in p.items():
        exponent, scalar = 1, F(1)
        for r in reversed(word):
            scalar *= -exponent
            exponent -= r
        assert exponent == -3
        detector += coeff * scalar
    assert detector == F(1,4)
    print('SURVIVAL DETECTOR: negative Witt action gives P(t)=(1/4)t^(-3), nonzero without PBW independence')
    print('LIMIT: finite algebra checks; nonzero survival and ideal/Cartan disjointness require the stated presentation argument')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
