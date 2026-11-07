#!/usr/bin/env python3
"""Exact finite checks for the Sarizadeh manifold-question preregistration.

These check the source's Section 2 displayed pushforward/integral identity
and the finite-time and coefficient identities used by the proposed proof.
No finite calculation establishes minimality, historic behavior on a manifold,
Cantor-measure existence, or the homeomorphic-measures theorem.
"""
from fractions import Fraction as Q
from itertools import permutations, product
from math import factorial
import sympy as sp


def push(f, mass):
    out = [Q(0) for _ in f]
    for x, a in enumerate(mass):
        out[f[x]] += a
    return tuple(out)


def empirical(f, x, n):
    assert n > 0
    out = [Q(0) for _ in f]
    for _ in range(n):
        out[x] += Q(1, n)
        x = f[x]
    return tuple(out)


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), Q(0))


def source_duality():
    """The source's finite-time integral identity, all maps on <=4 states."""
    cases = 0
    for d in range(1, 5):
        masses = [tuple(Q(int(x == j)) for x in range(d)) for j in range(d)]
        masses += [tuple(Q(x+1, d*(d+1)//2) for x in range(d))]
        phi = tuple(Q((-1)**x*(x+1), x+2) for x in range(d))
        for f in product(range(d), repeat=d):
            for mass in masses:
                for n in range(1, 13):
                    avg = [Q(0) for _ in range(d)]
                    current = mass
                    for _ in range(n):
                        avg = [a+b/n for a,b in zip(avg,current)]
                        current = push(f,current)
                    lhs = dot(phi,avg)
                    rhs = sum((mass[x]*dot(phi,empirical(f,x,n)) for x in range(d)),Q(0))
                    assert lhs == rhs
                    cases += 1
    print(f'FIDELITY: source Section 2 pushforward/integral identity: {cases} exact cases')


def conjugacy_and_shift():
    cases = 0
    for d in range(1, 5):
        for f in permutations(range(d)):
            for h in permutations(range(d)):
                inv = tuple(h.index(y) for y in range(d))
                g = tuple(inv[f[h[x]]] for x in range(d))
                for x in range(d):
                    for n in range(1, 13):
                        ef = empirical(f,h[x],n)
                        assert empirical(g,x,n) == push(inv,ef)
                        fnx = h[x]
                        for _ in range(n): fnx = f[fnx]
                        a = empirical(f,f[h[x]],n)
                        for y in range(d):
                            assert a[y]-ef[y] == Q(int(y==fnx)-int(y==h[x]),n)
                        cases += 1
    print(f'ROUTE: empirical conjugacy and endpoint shift: {cases} exact cases')


def exact_geometric_mixture():
    C = sp.symbols('C', positive=True)
    # Denominator C-1 is nonzero under the required C>1.
    a = (C-1)/(C+1)
    assert sp.cancel(a*(1+2/(C-1))-1) == 0
    checks = 0
    for den in range(1, 9):
        for num in range(den+1, 2*den+5):
            c = Q(num,den)
            w = lambda k: (c-1)/(c+1)*c**(-abs(k))
            for r in range(0, 13):
                # Infinite tail is summed algebraically, never truncated away.
                tail = 2*c**(-r)/(c+1)
                assert sum((w(k) for k in range(-r,r+1)),Q(0))+tail == 1
            for k in range(-24,25):
                assert w(k)>0
                assert Q(1,c) <= w(k-1)/w(k) <= c
                checks += 1
            # Exact wrapped infinite weights on each finite cyclic orbit.
            # Class j contains all integer k congruent to j modulo d.
            for d in range(1, 8):
                base = (c-1)/(c+1)
                q = Q(1,c)
                cycle = [base*(1+q**d)/(1-q**d)]
                cycle += [base*(q**j+q**(d-j))/(1-q**d) for j in range(1,d)]
                assert sum(cycle,Q(0)) == 1
                for j in range(d):
                    assert Q(1,c)*cycle[j] <= cycle[(j-1)%d] <= c*cycle[j]
                    checks += 1
    print(f'ROUTE: geometric normalization, exact tails, shifted weights: {checks} exact cases')


def oscillation_model():
    # A two-sided symbolic-shift model, not a torus homeomorphism example.
    # Each new constant block dominates the entire previous prefix.
    n = ones = 0
    for stage in range(1, 25):
        length = (stage+1)*(n+1)
        bit = stage % 2
        ones += bit*length
        n += length
        avg = Q(ones,n)
        if bit: assert avg > 1-Q(1,stage+1)
        else: assert avg < Q(1,stage+1)
    print('ROUTE: alternating symbolic blocks: 24 exact endpoint inequalities')


if __name__ == '__main__':
    source_duality()
    conjugacy_and_shift()
    exact_geometric_mixture()
    oscillation_model()
    print('LIMIT: finite checks support identities only; manifold existence remains unformalized')
    print('ALL_CHECKS_OK')
