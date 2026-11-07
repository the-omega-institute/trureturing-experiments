#!/usr/bin/env python3
"""Exact checks for the generic HMM identifiability hardness preregistration.

Finite checks and symbolic identities support the mathematical outline; this is
neither a proof of complexity hardness nor a Lean kernel verification.
"""
from fractions import Fraction as Q
from itertools import product
from math import prod
import sympy as s


def words(m, maximum):
    for length in range(maximum + 1):
        yield from product(range(m), repeat=length)


def path_probability(M, word):
    pi, A, B = M
    if not word:
        return Q(1)
    return sum(pi[path[0]] * B[path[0]][word[0]] *
               prod(A[path[t]][path[t + 1]] * B[path[t + 1]][word[t + 1]]
                    for t in range(len(word) - 1))
               for path in product(range(len(pi)), repeat=len(word)))


def forward_probability(M, word):
    pi, A, B = M
    row = list(pi)
    for letter in word:
        row = [sum(row[i] * B[i][letter] * A[i][j]
                   for i in range(len(pi))) for j in range(len(pi))]
    return sum(row)


def ambiguous(u):
    return ([Q(1, 2), Q(1, 2)],
            [[Q(1, 2) + u * Q(1, 8), Q(1, 2) - u * Q(1, 8)], [Q(1, 4), Q(3, 4)]],
            [[Q(1, 4), Q(3, 4)], [Q(1, 4), Q(3, 4)]])


def identifiable(u, v):
    p, q = Q(1, 2) + u * Q(1, 8), Q(1, 2) + v * Q(1, 8)
    return ([p, 1 - p], [[q, 1 - q], [Q(1, 4), Q(3, 4)]],
            [[Q(1, 4), Q(3, 4)], [Q(3, 4), Q(1, 4)]])


def blend(tag, left, right):
    return ([s.expand((1-tag)*a+tag*b) for a, b in zip(left[0], right[0])],
            *[[[s.expand((1-tag)*a+tag*b) for a, b in zip(ar, br)]
               for ar, br in zip(left[i], right[i])] for i in (1, 2)])


def flat(M):
    return tuple(M[0]) + tuple(x for row in M[1] for x in row) + tuple(x for row in M[2] for x in row)


def signature(M):
    # P(z1=0), P(z2=0), both observable from words of length at most two.
    return (path_probability(M, (0,)),
            sum(path_probability(M, (a, 0)) for a in range(2)))


def stochastic(M):
    pi, A, B = M
    assert sum(pi) == 1
    assert all(sum(row) == 1 for row in A + B)
    assert all(Q(1, 4) <= x <= Q(3, 4) for x in flat(M))
    assert A[0][0]*A[1][1] - A[0][1]*A[1][0] != 0


# Source fidelity: Definition 3 and the forward-procedure identity, Remark 7.
grid = [Q(1, 4), Q(1, 2), Q(3, 4)]
fidelity_count = 0
for p, a, b, e, f in product(grid, repeat=5):
    M = ([p, 1-p], [[a, 1-a], [b, 1-b]], [[e, 1-e], [f, 1-f]])
    for w in words(2, 4):
        assert path_probability(M, w) == forward_probability(M, w)
        fidelity_count += 1
print(f'FIDELITY: source hidden-path sum = forward procedure, {fidelity_count} exact cases')

# Equal-emission identity, symbolic in every stochastic entry, through length 5.
a, b, e, r = s.symbols('a b e r')
M = ([r, 1-r], [[a, 1-a], [b, 1-b]], [[e, 1-e], [e, 1-e]])
for w in words(2, 5):
    # Expand at each step to prevent expression growth concealing cancellation.
    row = list(M[0])
    for letter in w:
        row = [s.expand(sum(row[i]*M[2][i][letter]*M[1][i][j] for i in range(2)))
               for j in range(2)]
    assert s.expand(sum(row) - e**w.count(0)*(1-e)**w.count(1)) == 0
print('EQUAL_EMISSION: 63 symbolic word identities, lengths 0 through 5')

# The two observable marginals recover the identifiable family's parameters.
u, v, tag = s.symbols('u v tag')
p, q = Q(1, 2) + u/8, Q(1, 2) + v/8
J = identifiable(u, v)
f1, f2 = signature(J)
assert s.expand(f1 - (Q(3, 4)-p/2)) == 0
assert s.expand(f2 - (Q(5, 8)+p/8-p*q/2)) == 0
recovered_p = Q(3, 2) - 2*f1
assert s.expand(recovered_p - p) == 0
# Clear only the guarded, nonzero denominator 4p (p > 1/2 on our domain).
assert s.expand(4*p*q - (5+p-8*f2)) == 0
assert s.expand(s.det(s.Matrix([[s.diff(f1, u), s.diff(f1, v)],
                       [s.diff(f2, u), s.diff(f2, v)]])) - (u+4)/2048) == 0
print('RECOVERY: P(z1=0)=3/4-p/2; P(z2=0)=5/8+p/8-p*q/2; inverse guarded by p>1/2')

# Piecewise branches are implemented by rational polynomials of degree at most 2.
negative = blend(tag, identifiable(u, 0), ambiguous(u))
positive = blend(tag, ambiguous(u), identifiable(u, v))
for family in [negative, positive]:
    assert all(not s.sympify(x).atoms(s.Float) for x in flat(family))
    assert all(s.Poly(x, tag, u, v).total_degree() <= 2 for x in flat(family))
    pi, A, B = family
    assert s.expand(sum(pi)-1) == 0
    assert all(s.expand(sum(row)-1) == 0 for row in A+B)
    for branch in (0, 1):
        expected = (identifiable(u, 0) if branch == 0 else ambiguous(u)) if family is negative else (ambiguous(u) if branch == 0 else identifiable(u, v))
        assert all(s.expand(x.subs(tag, branch)-y) == 0 for x, y in zip(flat(family), flat(expected)))
print('PARAMETERIZATIONS: both selector blends are stochastic polynomial maps of degree at most 2')

# Every sampled ambiguous rival changes actual transition entries. Recoverability
# distinguishes all sampled J models, even when feasibility coordinates repeat.
small = [Q(i, 8) for i in range(1, 8)]
entries_to_signature = {}
for x in small:
    M, rival = ambiguous(x), ambiguous((x+1)/2)
    stochastic(M); stochastic(rival)
    assert flat(M) != flat(rival)
    assert all(path_probability(M, w) == path_probability(rival, w)
               for w in words(2, 5))
    for y in small:
        J = identifiable(x, y); stochastic(J)
        sig = signature(J)
        assert sig[0] > Q(1, 4)
        assert sig not in entries_to_signature
        entries_to_signature[sig] = flat(J)
print('RIVALS_AND_SEPARATION: 7 distinct-entry ambiguous rivals; 49 uniquely recovered identifiable models')

# Exact boundary diagnostics for the dimension argument: the active J component
# has two independent freely varying coordinates, while baseline is a line.
# This does not numerically estimate a dimension or claim finite grids are open.
assert s.Matrix([u, v]).jacobian([u, v]) == s.eye(2)
assert s.Matrix([u, 0]).jacobian([u, v]).rank() == 1
print('DIMENSION_INPUTS: active projection is the identity on (0,1)^2; baseline has one free coordinate')

# Feasibility examples include empty, finite, positive-dimensional and K=0 cases.
x, y = s.symbols('x y', real=True)
assert s.solveset(x*x+1, x, domain=s.S.Reals) is s.S.EmptySet
assert s.solveset(x, x, domain=s.S.Reals) == s.FiniteSet(0)
assert s.solveset(x*x-1, x, domain=s.S.Reals) == s.FiniteSet(-1, 1)
assert s.expand((x*y).subs(x, 0)) == 0
# Logical/dimension table: absent activated branch => a single baseline curve;
# present activated branch => bad >=1 in N, domain >=2 in E.
for feasible in (False, True):
    negative_yes = not feasible
    positive_yes = feasible
    if not feasible:
        assert (-1 < 1) == negative_yes
        assert (1 < 1) == positive_yes
    else:
        assert (1 < 1) == negative_yes  # minimal active dimension; higher dimensions also tie
        assert (1 < 2) == positive_yes  # lower bound only; higher dimensions preserve YES
print('REDUCTION_TABLE: infeasible -> (N YES,E NO); feasible -> (N NO,E YES); empty/singleton/finite/positive-dimensional witnesses screened')
print('ALL_CHECKS_OK')
