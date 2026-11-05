#!/usr/bin/env python3
"""Exact finite checks for Campolina--Mailybaev Remark 2; not an unbounded proof."""
from collections import defaultdict
from fractions import Fraction
from math import gcd
import sys
import sympy as sp

X = sp.Symbol('X')
BOUND = 60


def poly(n, m):
    return sp.Poly(X**n - X**m - 1, X, domain=sp.ZZ)


def signature(p):
    return tuple(map(int, p.all_coeffs()))


def table(a, b):
    return [((1, b), (-1, a)), ((-1, a), (1, b)),
            ((1, b-a), (-1, -a)), ((-1, -a), (1, b-a)),
            ((1, -b), (1, a-b)), ((1, a-b), (1, -b))]


def fidelity():
    # Source Table 1(a), and all twelve ordered triads in Table 1(c).
    dyadic = [(Fraction(2), Fraction(-1)), (Fraction(-1), Fraction(2)),
              (Fraction(1, 2), Fraction(1, 2))]
    assert all(p + q == 1 for p, q in dyadic)
    plastic = poly(3, 1)
    assert sp.expand((X**3-X-1)*(X**2-X+1)) == X**5-X**4-1
    rows = table(1, 3) + table(4, 5)
    assert len(set(rows)) == 12
    for ((s, i), (t, j)) in rows:
        k = max(0, -i, -j)
        assert sp.Poly(s*X**(i+k)+t*X**(j+k)-X**k, X).rem(plastic).is_zero
    # Exact cutoff calculation used by the source; no approximate root locations.
    x = Fraction(21, 20)
    assert x**63 - x**62 - 1 > 0
    print('FIDELITY: dyadic 3 triads; plastic 12 ordered triads; p_(62,63)(21/20)>0')


def main():
    print('Python:', sys.version.split()[0], 'SymPy:', sp.__version__)
    fidelity()
    groups = defaultdict(list)
    reducible = 0
    count = 0
    zero_middle = 0
    for n in range(1, BOUND+1):
        for m in range(n):
            f = poly(n, m)
            content, factors = sp.factor_list(f)
            assert content == 1
            reconstructed = sp.Poly(content, X, domain=sp.ZZ)
            cyclo = sp.Poly(1, X, domain=sp.ZZ)
            noncyclo = []
            for h, e in factors:
                assert h.is_irreducible
                reconstructed *= h**e
                if h.is_cyclotomic:
                    cyclo *= h**e
                else:
                    noncyclo.append((h, e))
            assert reconstructed == f
            assert len(noncyclo) == 1 and noncyclo[0][1] == 1
            q = noncyclo[0][0]
            assert q.LC() == 1
            # The unique real root >1 belongs to q: f(1)<0, f(2)>=0,
            # and cyclo has no positive real root >1. Count roots exactly.
            assert q.count_roots(1, sp.oo) == 1
            if m == 0:
                zero_middle += 1
                assert cyclo == sp.Poly(1, X) and q == f
            else:
                g = gcd(n, m)
                exceptional_residue = (n//g % 6, m//g % 6) in [(1, 2), (5, 4)]
                expected_c = sp.Poly(X**(2*g)-X**g+1, X) if exceptional_residue else sp.Poly(1, X)
                assert cyclo == expected_c
                if exceptional_residue:
                    reducible += 1
                    assert q.nth(0) == -1
                    assert all(q.nth(i) == 0 for i in range(1, g))
                    assert q.nth(g) == -1
            groups[signature(q)].append((n, m))
            count += 1
        if n % 10 == 0:
            print('FACTORED_THROUGH:', n, flush=True)
    shared = []
    for pairs in groups.values():
        if len(pairs) > 1:
            assert len(pairs) == 2
            first, second = sorted(pairs)
            g = gcd(*first)
            assert first == (3*g, g) and second == (5*g, 4*g)
            shared.append(g)
    assert sorted(shared) == list(range(1, BOUND//5+1))
    print('EXACT_FACTORIZATION:', count, 'polynomials; 1<=n<=60; 0<=m<n;',
          zero_middle, 'm=0 Eisenstein cases;', reducible, 'nontrivial cyclotomic parts')
    print('COMMON_POSITIVE_ROOT_PAIRS:', len(shared), 'exactly (3g,g),(5g,4g), g=1..12')
    # Independently verify both polynomial identities for every exceptional g in range.
    for g in shared:
        assert poly(3*g, g) * sp.Poly(X**(2*g)-X**g+1, X) == poly(5*g, 4*g)
    # Combinatorial sign/exponent graph of one relation. Its projected generators
    # are a,b; each hyperedge connects both signs, so residue classes mod gcd(a,b)
    # are the components. Check the finite cyclic analogue independently by BFS.
    trials = 0
    for b in range(1, BOUND+1):
        for a in range(b):
            for period in [b, 2*b, 2*b+1]:
                vertices = [(s, i) for s in [-1, 1] for i in range(period)]
                adjacency = {v:set() for v in vertices}
                for s, i in vertices:
                    triples = [[(s, (i+b)%period), (s, (i+a)%period), (s, i)],
                               [(s, (i+b)%period), (-s, (i+a)%period), (s, i)],
                               [(s, (i+b)%period), (s, (i+a)%period), (-s, i)]]
                    for triple in triples:
                        for u in triple:
                            for v in triple:
                                adjacency[u].add(v)
                todo=[(1,0)]; seen=set(todo)
                while todo:
                    for v in adjacency[todo.pop()]:
                        if v not in seen: seen.add(v);todo.append(v)
                d=gcd(gcd(a,b),period)
                expected={(s,i) for s,i in vertices if i%d==0}
                assert seen == expected
                trials += 1
    print('FINITE_CYCLIC_COMPONENT_CHECKS:', trials)
    print('SCOPE: finite exact evidence only; the unrestricted assertion is not proved here')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
