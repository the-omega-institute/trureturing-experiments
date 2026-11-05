#!/usr/bin/env python3
"""Exact checks for Goldberg arXiv:1407.7430v3, Elphick Conjecture 6.2.

Uses Python integers and SymPy exact expressions only. No numerical eigensolver.
Run with Python 3 and SymPy. This is experimental evidence, not a Lean proof.
"""
import sys
import sympy as S

x = S.Symbol('x')


def adjacency(n, edges):
    A = S.zeros(n)
    edges = list(edges)
    assert len(set(tuple(sorted(e)) for e in edges)) == len(edges)
    for a, b in edges:
        assert 0 <= a < n and 0 <= b < n and a != b
        A[a, b] = A[b, a] = 1
    seen, todo = {0}, [0]
    while todo:
        a = todo.pop()
        for b in range(n):
            if A[a, b] and b not in seen:
                seen.add(b)
                todo.append(b)
    assert len(seen) == n
    assert A == A.T and S.trace(A) == 0
    assert S.trace(A * A) == 2 * len(edges)
    return A


def tree(t, extra=False):
    edges = [(0, i) for i in range(1, t + 1)]
    edges += [(i, t + i) for i in range(1, t + 1)]
    if extra:
        edges += [(0, 2 * t + 1)]
    A = adjacency(2 * t + 1 + int(extra), edges)
    assert len(edges) == A.rows - 1  # connected, hence a tree
    return A, len(edges)


def same_poly(A, expected):
    assert S.Poly(A.charpoly(x).as_expr() - expected, x).is_zero


# Fidelity: Section 1 trace identity and Section 3 complete-bipartite
# energy E=2 sqrt(pq), equality in Corollary 3.2 E >= 2m/lambda_1.
# The source's adjacent text says m=2pq; actual unordered edge count is pq.
for p in range(1, 5):
    for q in range(1, 6):
        A = adjacency(p + q, ((a, p + b) for a in range(p) for b in range(q)))
        same_poly(A, x ** (p + q - 2) * (x*x - p*q))
        rho = S.sqrt(p*q)
        E = 2*rho
        assert S.simplify(E - 2*p*q/rho) == 0
print('FIDELITY: 20 complete bipartite graphs; exact charpoly, trace(A^2)=2m,')
print('E=2*sqrt(p*q)=2m/rho. Source typo m=2pq is excluded; actual m=pq.')

# Actual adjacency characteristic polynomials, retaining all multiplicities.
for t in list(range(1, 13)) + [15, 16, 20]:
    A, m = tree(t)
    same_poly(A, x*(x*x-t-1)*(x*x-1)**(t-1))
print('T_t: exact full adjacency charpolys at t=1..12,15,16,20 agree.')

# Independent rational eigenbasis of the actual 31-by-31 matrix, not just
# a supplied polynomial: 28 arm differences, one zero mode, two hub modes.
A, m = tree(15)
vectors, eigenvalues = [], []
for sign in [-1, 1]:
    for i in range(1, 15):
        v = S.zeros(31, 1)
        v[i], v[15], v[15+i], v[30] = 1, -1, sign, -sign
        vectors.append(v)
        eigenvalues.append(sign)
v = S.zeros(31, 1)
v[0] = 1
for i in range(1, 16):
    v[15+i] = -1
vectors.append(v)
eigenvalues.append(0)
for lam in [-4, 4]:
    v = S.ones(31, 1)
    v[0] = 15
    for i in range(1, 16):
        v[i] = lam
    vectors.append(v)
    eigenvalues.append(lam)
P = S.Matrix.hstack(*vectors)
assert P.shape == (31, 31) and P.rank() == 31
assert A*P == P*S.diag(*eigenvalues)
rho, E = max(eigenvalues), sum(abs(a) for a in eigenvalues)
bound = S.Integer(2*m)/S.sqrt(rho)
assert (A.rows, m, rho, E, bound) == (31, 30, 4, 36, 30)
assert E > bound
print('WITNESS: n=31 m=30 rho=4 E=36 proposed_bound=30; basis rank=31.')

# Symbolic invariant-layer quotient, arithmetic gap, and unbounded ratio.
t, k = S.symbols('t k', integer=True, positive=True)
Q = S.Matrix([[0,t,0], [1,0,1], [0,1,0]])
assert S.expand(Q.charpoly(x).as_expr()-x*(x*x-t-1)) == 0
gap = 2*(k**4+k**2-2)-4*(k**4-1)/k
assert S.simplify(k*gap/2-(k-1)*(k**3*(k-1)-2)) == 0
ratio = k*(k**4+k**2-2)/(2*(k**4-1))
assert S.simplify(ratio-k/2-k/(2*(k*k+1))) == 0
for kval in range(2, 101):
    assert gap.subs(k,kval) > 0
    assert ratio.subs(k,kval) > S.Rational(kval,2)
print('FAMILY: exact gap factorization; ratio=k/2+k/(2*(k^2+1))>k/2.')
print('Checked strict gap for k=2..100; symbolic formulas cover the proposed route.')

# Robustness: U_t adds a central leaf and has no zero eigenvalue.
for tt in list(range(1, 9)) + [16, 20]:
    U, mm = tree(tt, True)
    same_poly(U, (x*x-1)**(tt-1)*(x**4-(tt+2)*x*x+1))
    assert U.det() != 0
r = (S.sqrt(t+4)+S.sqrt(t))/2
s = (S.sqrt(t+4)-S.sqrt(t))/2
assert S.simplify(r*s-1) == 0
assert S.simplify(r+s-S.sqrt(t+4)) == 0
for tt in [16, 17, 20, 100, 10000]:
    rr = r.subs(t,tt)
    EE = 2*(tt-1+S.sqrt(tt+4))
    assert S.simplify(rr-4).is_nonnegative
    assert S.simplify(EE-(2*tt+2)).is_nonnegative
    assert S.simplify(EE-(2*tt+1)).is_positive  # 2m/sqrt(rho) <= 2t+1 < E
print('U_t: exact charpolys and nonzero determinants at t=1..8,16,20;')
print('positive reciprocal exceptional roots; radical bounds checked at t=16,17,20,100,10000.')
print('ALL_CHECKS_OK')
