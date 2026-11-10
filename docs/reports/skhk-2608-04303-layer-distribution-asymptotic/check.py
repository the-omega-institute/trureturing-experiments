#!/usr/bin/env python3
"""Exact finite fidelity and identity checks; no asymptotic certification."""
from fractions import Fraction as Q
from itertools import product
from math import factorial
import sympy as sp


def layers(parents):
    n = len(parents) + 1
    adj = [set() for _ in range(n)]
    for child, parent in enumerate(parents, 1):
        adj[child].add(parent)
        adj[parent].add(child)
    alive, result, k = set(range(n)), [-1] * n, 0
    while alive:
        leaves = {v for v in alive if len(adj[v] & alive) <= 1}
        assert leaves, 'a nonempty finite tree has a leaf or is a singleton'
        for v in leaves:
            result[v] = k
        alive -= leaves
        k += 1
    return result


def source_l1(n):
    return (-Q(1, 2) + Q(1, n) - Q(1, n - 1) + Q(2, factorial(n))
            + sum((Q(1, factorial(v)) for v in range(2, n)), Q(0))
            - Q(1, n * (n - 1)) * sum((Q(1, factorial(v)) for v in range(n - 2)), Q(0))
            + Q(2, factorial(n)) * sum((Q(factorial(j), factorial(v))
                for j in range(1, n - 1) for v in range(j)), Q(0)))


print('Exact arithmetic only: integers, Fraction, symbolic exp(1).')
for n in range(3, 9):
    total, c0, c1 = 0, 0, 0
    for parents in product(*(range(i) for i in range(1, n))):
        ls = layers(parents)
        total += 1
        c0 += ls.count(0)
        c1 += ls.count(1)
    assert total == factorial(n - 1)
    e0, e1 = Q(c0, n * total), Q(c1, n * total)
    assert e0 == Q(1, 2) + Q(1, n * (n - 1))
    assert e1 == source_l1(n)
    print(f'Source Theorem 7: n={n}, histories={total}, ell0(n)={e0}, ell1(n)={e1}: OK')

# z=exp(-y): B_k=integral_0^1 S_k(z) dz.
z = sp.Symbol('z', real=True)
S0, S1, S2 = sp.Integer(1), 1-z, 1-z*sp.exp(1-z)
Bs = [sp.integrate(f, (z, 0, 1)) for f in (S0, S1, S2)]
assert [sp.simplify(a-b) for a,b in zip(Bs,[1, sp.Rational(1,2), 3-sp.E])] == [0,0,0]
print('Source Theorem 2 / small-k evaluation: B0=1, B1=1/2, B2=3-exp(1): OK')
for k, (a,b) in enumerate([(S0,S1),(S1,S2)]):
    mass = sp.simplify(Bs[k] - Bs[k+1])
    positive_integral = sp.integrate(a*b, (z, 0, 1))
    assert sp.simplify(mass-positive_integral) == 0
    print(f'Positive integral identity: k={k}, ell_k={mass}: OK')

# Exact formal Taylor coefficients of the source recurrence at y=0.
# Coefficients of exp(a) satisfy n*b_n=sum(j*a_j*b_(n-j)).
N = 18

def mul(a,b):
    return [sum((a[j]*b[d-j] for j in range(d+1)),Q(0)) for d in range(N+1)]

def next_s(a):
    exponent = [Q(0)] + [-a[j-1]/j for j in range(1,N+1)]
    e = [Q(1)] + [Q(0)]*N
    for d in range(1,N+1):
        e[d] = sum((j*exponent[j]*e[d-j] for j in range(1,d+1)),Q(0))/d
    return [Q(0)] + [-x for x in e[1:]]

s = [[Q(1)] + [Q(0)]*N]
for _ in range(9):
    s.append(next_s(s[-1]))
for k in range(9):
    a,b = s[k],s[k+1]
    ab = mul(a,b)
    derivative = [(j+1)*b[j+1] for j in range(N)]
    # s'_(k+1)=s_k*(1-s_(k+1)), and the weighted integration-by-parts residual.
    assert derivative == [a[j]-ab[j] for j in range(N)]
    assert [a[j]-b[j]-ab[j] for j in range(N)] == [derivative[j]-b[j] for j in range(N)]
    print(f'Exact Taylor identity residual: k={k}, coefficients y^0..y^{N-1}: zero')
print('Scope: finite-n fidelity; exact k=0,1 integrals; finite Taylor coefficients, not a proof of the asymptotic or of the fixed-k limit.')
print('ALL_CHECKS_OK')
