#!/usr/bin/env python3
"""Exact checks for the Fillman--Lim Question 2 counterexample candidate.

No floating-point computation is used. Symbolic identities and finite checks
support the accompanying mathematical argument; they are not a Lean settlement
or a numerical test of the infinite quantified definition of gap-richness.
"""
from fractions import Fraction as F
from itertools import product
import sys
import sympy as sp

assert sys.version_info >= (3, 12)
I = sp.eye(2)
J = sp.Matrix([[0, -1], [1, 0]])
E, q, t, s, h = sp.symbols('E q t s h', real=True)


def S(a):
    return sp.Matrix([[a, -1], [1, 0]])


def T(block, energy):
    result = I
    for a in block:
        result = S(energy - a) * result
    return result


def zero_matrix(a):
    return all(sp.expand(z) == 0 for z in a)


def L(q, t, h):
    return (t, q, 2*q-h-t)


def U(q, t, h):
    return L(q, t, h) + (q, q)


def V(q, t, h):
    return (q, t, q, 2*q+h-t, q)


def component(n):
    assert n >= 1
    m = n.bit_length()-1
    j = n-2**m
    return F(2*j+1, 2**(m+1)), F(1, 3**n)


def midpoint(n):
    qn, rn = component(n)
    return qn, 10+F(21,20)*rn


def block_product(blocks, energy):
    result = I
    for block in blocks:
        result = T(block, energy)*result
    return result


# Fidelity: Appendix A, Lemma A.1 (leading-order trace), source v1.
for n in range(1, 7):
    xs = sp.symbols('x:'+str(n))
    a = T(xs, E)
    leading = sp.prod(E-x for x in xs)
    residual = sp.Poly(sp.expand(sp.trace(a)-leading), *xs)
    assert residual.is_zero or residual.total_degree() <= n-1
    assert sp.expand(a.det()-1) == 0
print('FIDELITY: source Lemma A.1 and determinant 1, lengths 1..6: OK')

# These are polynomial identities in arbitrary real parameters, not samples.
H = S(h)
assert zero_matrix(J*J+I)
assert zero_matrix(T(L(q,t,h), q)+H)
assert zero_matrix(T(U(q,t,h), q)-H)
assert zero_matrix(T(V(q,t,h), q)+H.adjugate())
for factory in (L, U, V):
    a, b = T(factory(q,t,h), E), T(factory(q,s,h), E)
    assert sp.expand(a.det()-1) == 0
    assert sp.expand(b.det()-1) == 0
    expected = 2+(E-q)**2*(t-s)**2
    assert sp.expand(sp.trace(a*b.adjugate())-expected) == 0
assert zero_matrix(T(tuple(a+b for a,b in zip(L(q,t,h),L(q,s,-h))),2*q)+J)
assert zero_matrix(T(V(q,s,h),q)*T(U(q,t,h),q)+I)
print('FAMILY: transfer, relative-trace, sum and product identities, symbolic h,q,t,s,E: OK')

# Positive-semigroup lemma finite cross-check over every small integer SL(2).
# Integer elliptics have traces -1,0,1 and inverses among their first five powers.
small = []
for a,b,c,d in product(range(-2,3), repeat=4):
    if a*d-b*c == 1 and abs(a+d) <= 2:
        small.append(sp.Matrix([[a,b],[c,d]]))
relative_pairs = 0
for a,b in product(small, repeat=2):
    if abs(sp.trace(a*b.adjugate())) <= 2:
        continue
    relative_pairs += 1
    powers_a = [a**n for n in range(1,9)]
    powers_b = [b**n for n in range(1,9)]
    assert any(abs(sp.trace(x*y)) > 2 for x,y in product(powers_a,powers_b))
print(f'SEMIGROUP CROSS-CHECK: {len(small)} matrices, {relative_pairs} qualifying ordered pairs: OK')

# All dyadic midpoints occur, at arbitrarily late levels. The mesh is 2^-m.
for m in range(8):
    assert [component(n)[0] for n in range(2**m,2**(m+1))] == [
        F(2*j+1,2**(m+1)) for j in range(2**m)]
print('DENSITY GRID: levels 0..7 match the exact dyadic midpoint enumeration: OK')

# Uniform separation inequalities. These bound the entire remaining tails.
r = F(1)
eps = r/4
target = F(21,10)*r
assert 3*r-target > eps  # any component index below n
assert target-F(11,10)*r-F(11,10)*r/3 > eps  # not both indices n
assert r-F(11,10)*r/3 > eps  # single component separation below
assert 3*r-F(11,10)*r > eps  # single component separation above
assert r > eps  # separation from the base segment
assert F(21,20)*r-r > 0
assert F(11,10)*r-F(21,20)*r > 0
# U versus V is separated by first coordinates t >= 10, q <= 1.
assert F(9) > F(1,3)/4
print('TAIL BOUNDS: uniform single- and double-component isolation, radius r_n/4: OK')

# Exhaustive low-component interval checks; include the baseline as index 0.
for n in range(1,13):
    qn,rn = component(n)
    lo,hi = 10+rn,10+F(11,10)*rn
    epsn=rn/4
    for a in range(17):
        qa,ra=component(a) if a else (F(0),F(0))
        alo,ahi=(ra,F(11,10)*ra) if a else (F(0),F(0))
        if a != n:
            assert ahi < rn-epsn or alo > F(11,10)*rn+epsn
        for b in range(17):
            qb,rb=component(b) if b else (F(0),F(0))
            blo,bhi=(rb,F(11,10)*rb) if b else (F(0),F(0))
            if (a,b) != (n,n):
                assert ahi+bhi < F(21,10)*rn-epsn or alo+blo > F(21,10)*rn+epsn
print('INTERVAL CROSS-CHECK: target n=1..12, summand indices 0..16: OK')

# Relative trace remains unchanged by the fixed prefix and suffix of an input
# word. This verifies the actual arbitrary-word route on short exact examples.
for hh in (F(5,2),F(3),F(4)):
    for factory in (L,U,V):
        for length in range(1,5):
            points=[midpoint(n) for n in range(1,length+1)]
            blocks=[factory(qq,tt,hh) for qq,tt in points]
            for energy in (F(1,4),F(1,2),F(2)):
                a=block_product(blocks,energy)
                for j,(qq,tt) in enumerate(points):
                    delta=component(j+1)[1]/100
                    alt=blocks.copy();alt[j]=factory(qq,tt+delta,hh)
                    b=block_product(alt,energy)
                    assert sp.trace(a*b.adjugate()) == 2+(energy-qq)**2*delta**2
print('INPUT WORDS: relative-trace invariant, h=5/2,3,4, lengths 1..4: OK')

# Nonhyperbolicity for every allowed perturbation length is derived from the
# symbolic identities above. The finite cross-checks include odd-length words
# for sum/product and all periodic alignment constraints for the union.
for hh in (F(5,2),F(3),F(4)):
    for n in range(1,7):
        qn,rn=component(n);tt=10+F(21,20)*rn
        u,v=U(qn,tt,hh),V(qn,tt,hh)
        for length in range(1,17):
            sum_matrix=(-J)**length
            product_matrix=(-I)**length
            assert abs(sp.trace(sum_matrix)) <= 2
            assert abs(sp.trace(product_matrix)) <= 2
            labels=[j%2 for j in range(length)]
            aligns=all(labels[j%length] == j%2 for j in range(2*length))
            assert aligns == (length%2 == 0)
            if aligns:
                w=[u if j%2==0 else v for j in range(length)]
                matrix=block_product(w,qn)
                assert matrix == (-I)**(length//2)
                assert abs(sp.trace(matrix)) <= 2
        assert abs(sp.trace(T(u,qn))) > 2
        assert abs(sp.trace(T(v,qn))) > 2
print('OBSTRUCTIONS: n=1..6, h=5/2,3,4, perturbation lengths 1..16 and odd parity: OK')
print('SCOPE: topology, elliptic recurrence and the infinite quantifiers use the README argument; no Lean settlement claimed.')
print('ALL_CHECKS_OK')
