#!/usr/bin/env python3
"""Exact finite checks for arXiv:2608.29574v1, Section 9 (first direction).

Run with Python >=3.12 and SymPy 1.14.0. No floating point is used.
These bounded checks are evidence for a preregistration, not a universal proof.
"""
from itertools import combinations_with_replacement
from random import Random
import sys
import sympy as s

assert sys.version_info >= (3, 12)
z, x = s.symbols('z x')
I = s.I
SEED = 260829574
rng = Random(SEED)
counts = {'diagonal': 0, 'jordan': 0, 'conjugated': 0,
          'cyclic': 0, 'inverse': 0, 'nilpotent_spaces': 0,
          'nonnilpotent_spaces': 0}


def clean(M):
    return M.applyfunc(s.expand)


def zero(M):
    return all(s.expand(v) == 0 for v in M)


def delta(A, T):
    return clean(A*T - T*A)


def poly(roots):
    return s.Poly(s.expand((2+I)*s.prod(z-r for r in roots)), z)


def action(P, A, T):
    out = s.zeros(A.rows)
    power = T
    for k in range(P.degree()+1):
        out += P.nth(k)*power
        power = delta(A, power)
    return clean(out)


def op(A):
    # Row-major vectorization: vec_row(AT-TA) = (A⊗I-I⊗Aᵀ)vec_row(T).
    n = A.rows
    return clean(s.kronecker_product(A, s.eye(n)) -
                 s.kronecker_product(s.eye(n), A.T))


def eval_op(P, L):
    out = s.zeros(L.rows)
    for c in P.all_coeffs():
        out = clean(out*L+c*s.eye(L.rows))
    return out


def word(roots, bound):
    distinct = list(dict.fromkeys(roots))
    for k in range(1, bound+1):
        for w in combinations_with_replacement(distinct, k):
            if s.expand(sum(w)) == 0:
                return w
    return None


def has_cycle(eigenvalues, P):
    n = len(eigenvalues)
    edges = [[s.expand(P.eval(eigenvalues[i]-eigenvalues[j])) == 0
              for i in range(n)] for j in range(n)]
    def visit(v, visited):
        for w in range(n):
            if edges[v][w] and (w in visited or visit(w, visited | {w})):
                return True
        return False
    return any(visit(v, {v}) for v in range(n))


def check_space(A, eigenvalues, P, roots, kind):
    n = A.rows
    L = op(A)
    Q = eval_op(P, L)
    basis = Q.nullspace()
    # Independently verify vectorization and every returned exact kernel vector.
    generic = s.zeros(n)
    coeffs = s.symbols(f'u0:{len(basis)}')
    for u, v in zip(coeffs, basis):
        B = s.Matrix(n, n, list(v))
        assert zero(action(P, A, B))
        assert zero(s.Matrix(n, n, list(L*v)) - delta(A, B))
        generic += u*B
    # All elements of a complex linear space are nilpotent iff the characteristic
    # polynomial of its generic element is identically x^n. A nonzero coefficient
    # polynomial takes a nonzero value over the infinite field ℂ.
    cp = s.Poly(generic.charpoly(x).as_expr()-x**n, x)
    all_nilpotent = all(s.expand(c) == 0 for c in cp.all_coeffs())
    if all_nilpotent:
        assert zero(generic**n)
        counts['nilpotent_spaces'] += 1
    else:
        assert word(roots, n) is not None
        counts['nonnilpotent_spaces'] += 1
    if kind in ['diagonal', 'conjugated']:
        # A fixed A can fail to realize a root word; compare its kernel with its
        # own support graph, and only use the global word condition as necessity.
        assert (not all_nilpotent) == has_cycle(eigenvalues, P)
    counts[kind] += 1
    return len(basis), all_nilpotent


def construction(roots, n):
    w = word(roots, n)
    if w is None:
        return
    k = len(w)
    a = [s.Integer(0)]
    for r in w[:-1]:
        a.append(s.expand(a[-1]+r))
    A = s.diag(*(a+[0]*(n-k)))
    T = s.zeros(n)
    for j in range(k):
        T[(j+1) % k, j] = 1
    P = poly(roots)
    assert zero(action(P, A, T))
    projector = s.diag(*([1]*k+[0]*(n-k)))
    assert T**k == projector
    assert T**(k*((n+k-1)//k)) == projector
    if all(r != 0 for r in w):
        S = T
        for h in range(5):
            assert zero(clean(S**k-s.prod(w)**h*projector))
            S = delta(A, S)
    counts['cyclic'] += 1


# Fidelity: Section 6, the source's displayed 3×3 example (roots 2 and -1).
A = s.diag(1, -1, 0)
T = s.Matrix([[0,1,0], [0,0,1], [1,0,0]])
S = delta(A,T)
assert S == s.Matrix([[0,2,0], [0,0,-1], [-1,0,0]])
assert delta(A,S) == s.Matrix([[0,4,0], [0,0,1], [1,0,0]])
assert zero(delta(A,S)-S-2*T)
assert T.det() == 1 and S.det() == 2
print('FIDELITY: source Section 6 3x3 resonance example reproduced exactly')

# A genuine third-order, non-collinear-root case, with sharp dimension threshold.
roots3 = [s.Integer(1), I, -1-I]
P3 = s.Poly(s.prod(z-r for r in roots3).expand(), z)
A3 = s.diag(0, 1, 1+I)
T3 = s.Matrix([[0,0,1], [1,0,0], [0,1,0]])
assert P3.as_expr() == z**3-I*z-1+I
assert zero(action(P3, A3, T3))
assert T3**3 == s.eye(3)
assert zero(delta(A3,T3)**3 - (1-I)*s.eye(3))
assert word(roots3,2) is None and len(word(roots3,3)) == 3
print('THIRD_ORDER: P=z^3-i*z-1+i, T^3=I, [A,T]^3=(1-i)I; shortest word=3')

adversarial = [
    [0], [1], [-1], [1,2], [1,-1], [2,-1], [3,-1], [4,-1],
    roots3, [1,I,-1,-I], [1,1+I,2-I], [1,1,-1,-1],
    [1,I,-1-I,-1-I], [0,0,1], [1,1,1],
]
root_pool = list(map(s.sympify, [-2,-1,0,1,2,I,-I,1+I,-1-I]))
random_roots = [[rng.choice(root_pool) for _ in range(rng.randint(1,4))]
                for _ in range(12)]
cases = [list(map(s.sympify,r)) for r in adversarial]+random_roots

for roots in cases:
    P = poly(roots)
    for n in range(1,5):
        eigenvalues = [rng.choice(root_pool) for _ in range(n)]
        check_space(s.diag(*eigenvalues), eigenvalues, P, roots, 'diagonal')
        construction(roots,n)
print(f'SMALL_CASES: {len(adversarial)} adversarial and {len(random_roots)} seeded random polynomials; n=1..4')

# The generic-kernel test is performed after a genuine integer similarity,
# independently of the entry-support computation for a diagonal A.
for n in range(1,5):
    for roots in [roots3, [0,0,1], [2,-1]]:
        eigenvalues = [s.Integer(j) for j in range(n)]
        U = s.eye(n)
        for i in range(n):
            for j in range(i+1,n):
                U[i,j] = rng.randint(-2,2)
        A = clean(U*s.diag(*eigenvalues)*U.inv())
        check_space(A,eigenvalues,poly(roots),roots,'conjugated')

# Non-diagonalizable A, several primary blocks, and repeated roots of P.
J2 = s.Matrix([[0,1],[0,0]])
J3 = s.Matrix([[0,1,0],[0,0,1],[0,0,0]])
jordan_cases = [
    (J2,[0,0]),
    (s.diag(J2,s.Matrix([[1]])),[0,0,1]),
    (s.diag(J2,J2+s.eye(2)),[0,0,1,1]),
    (s.diag(J3,s.Matrix([[2]])),[0,0,0,2]),
]
for A,eigenvalues in jordan_cases:
    for roots in [[0],[0,0],[1],[1,1],[1,-1],[1,1,-1,-1],roots3,[2,-1]]:
        check_space(A,eigenvalues,poly(list(map(s.sympify,roots))),list(map(s.sympify,roots)),'jordan')
L2 = op(J2)
assert len(L2.nullspace()) == 2
assert len((L2**2).nullspace()) == 3
print('REPEATED_ROOT_BOUNDARY: dim ker(ad_J2)=2, dim ker(ad_J2^2)=3; support criterion unchanged')

# Finite geometric inverse for nilpotent perturbations on Hom(V_mu,V_lambda).
for p in range(1,5):
    for q in range(1,5):
        Jp,Jq = s.zeros(p),s.zeros(q)
        for j in range(p-1): Jp[j,j+1] = 1
        for j in range(q-1): Jq[j,j+1] = 1
        D = s.kronecker_product(Jp,s.eye(q))-s.kronecker_product(s.eye(p),Jq.T)
        d = p+q-1
        E = s.eye(p*q)
        assert zero(D**d)
        # P(z)=z²-1, a=2, P(a)=3, B=(P(aI+D)-3I)/3.
        Q = (2*E+D)**2-E
        B = (Q-3*E)/3
        assert zero(B**d)
        inverse = sum(((-B)**h for h in range(d)),s.zeros(p*q))/3
        assert Q*inverse == E and inverse*Q == E
        counts['inverse'] += 1

# Scalar boundary and the sharp second-order arithmetic thresholds through n=4.
for roots in cases:
    P = poly(roots)
    assert (P.eval(0) == 0) == (word(roots,1) is not None)
for p in range(1,5):
    for q in range(1,5):
        if s.gcd(p,q) == 1:
            for n in range(1,5):
                assert (word([s.Integer(p),-s.Integer(q)],n) is not None) == (p+q <= n)
print('EXACT_COUNTS:', ', '.join(f'{k}={v}' for k,v in counts.items()))
print(f'ARITHMETIC: SymPy {s.__version__}; integers/rationals/Gaussian rationals; random seed={SEED}')
print('SCOPE: finite checks only; no universal proof or Lean proof is claimed')
print('ALL_CHECKS_OK')
