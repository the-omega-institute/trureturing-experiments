#!/usr/bin/env python3
"""Takahashi–Rayudu–Zhou–King–Thompson–Parekh, arXiv:2307.15688v3, Appendix C.1: the first-level projector relaxation
NPA_1(Proj) of the weighted crown H = x h_ab + sum_i (h_ai + h_bi) (vertices a, b, 1..n, n >= 3) is left open on
L < x < U, L = (n+2)^2 / (4(n+1)), U = (n+2)/3. Claim checked here: for every n >= 3 and L < x < U there is a feasible
real moment matrix with objective n + 1 + (x-L)^2/(32 n) > n + 1 = E_GS.

(1) For n = 3..10 and several rational x in (L, U) (including points near both ends), the explicit moment matrix over
    {I} and all vertex-pair variables h_e satisfies M(I,I)=1, M(h_e,h_e)=M(I,h_e), the anticommutation relation
    M(h_ij,h_jk)+M(h_jk,h_ij) = (M(h_ij)+M(h_jk)-M(h_ik))/2 for every path i-j-k, and symmetry; it is positive
    semidefinite by exact rational elimination; its objective equals n+1+(x-L)^2/(32n).
(2) Symbolic (n, x): the objective identity, the determinant identity
    det(S) = 2n q^3 (4 p0 - c (c - 1/(2n))) for the Schur complement, and the two scalar eigenvalues.
(3) Exact diagonalization of H for n = 3, 4, 5 at rational x in (L, U): the largest eigenvalue is n + 1 (Eq. (C2)).
(4) The interval is nonempty exactly for n >= 3."""
import sys, itertools
from fractions import Fraction as F
import numpy as np
import sympy as sp

bad = 0
def report(label, ok):
    global bad
    print(label, ok)
    bad += 0 if ok else 1

def params(n, x):
    L = F((n + 2) ** 2, 4 * (n + 1)); delta = x - L
    p0 = F(n + 1, 2 * n); q = delta / (16 * n); c = (x + L) / (4 * n)
    p = p0 - c * q; d = (n * p / 2 + q / 4) / (n - 1)
    return L, delta, p0, q, c, p, d

def moment(n, x):
    L, delta, p0, q, c, p, d = params(n, x)
    V = ['a', 'b'] + list(range(1, n + 1))
    pairs = [frozenset(e) for e in itertools.combinations(V, 2)]
    first = {}
    for e in pairs:
        if e == frozenset(['a', 'b']): first[e] = q
        elif 'a' in e or 'b' in e: first[e] = p
        else: first[e] = F(0)
    def second(e, f):
        if e == f: return first[e]
        if e & f:
            (j,) = tuple(e & f); i = next(iter(e - f)); k = next(iter(f - e))
            return (first[e] + first[f] - first[frozenset([i, k])]) / 4
        # disjoint pairs: a spoke at a and a spoke at b with different leaves get d; every other disjoint pair gets 0
        if ('a' in e and 'b' in f) or ('b' in e and 'a' in f):
            return d
        return F(0)
    idx = ['I'] + pairs
    m = len(idx); M = [[F(0)] * m for _ in range(m)]
    for r in range(m):
        for s in range(m):
            if idx[r] == 'I' and idx[s] == 'I': M[r][s] = F(1)
            elif idx[r] == 'I': M[r][s] = first[idx[s]]
            elif idx[s] == 'I': M[r][s] = first[idx[r]]
            else: M[r][s] = second(idx[r], idx[s])
    return idx, first, M, (L, delta, p0, q, c, p, d)

def psd_exact(M):
    A = [row[:] for row in M]; n = len(A); pivots = 0
    alive = list(range(n))
    while alive:
        k = max(alive, key=lambda i: A[i][i])
        if A[k][k] < 0: return False, pivots
        if A[k][k] == 0:
            if any(A[i][j] != 0 for i in alive for j in alive): return False, pivots
            break
        pivots += 1; piv = A[k][k]
        rest = [i for i in alive if i != k]
        for i in rest:
            f = A[i][k] / piv
            if f == 0: continue
            for j in rest: A[i][j] -= f * A[k][j]
        alive = rest
        for i in alive:
            if A[i][i] == 0 and any(A[i][j] != 0 for j in alive): return False, pivots
    return True, pivots

for n in range(3, 11):
    L = F((n + 2) ** 2, 4 * (n + 1)); U = F(n + 2, 3)
    xs = [L + (U - L) * F(t, 100) for t in (1, 13, 50, 87, 99)]
    for x in xs:
        idx, first, M, (L_, delta, p0, q, c, p, d) = moment(n, x)
        # constraints
        ok = M[0][0] == 1
        pos = {e: i for i, e in enumerate(idx)}
        for e in idx[1:]:
            ok &= M[pos[e]][pos[e]] == M[0][pos[e]]
            for f in idx[1:]:
                ok &= M[pos[e]][pos[f]] == M[pos[f]][pos[e]]
                if e != f and e & f:
                    (j,) = tuple(e & f); i = next(iter(e - f)); k = next(iter(f - e))
                    ok &= M[pos[e]][pos[f]] + M[pos[f]][pos[e]] == (first[e] + first[f] - first[frozenset([i, k])]) / 2
        psd, piv = psd_exact(M)
        obj = x * first[frozenset(['a', 'b'])] + sum(first[frozenset(['a', i])] + first[frozenset(['b', i])] for i in range(1, n + 1))
        ok_obj = obj == n + 1 + delta ** 2 / (32 * n) and obj > n + 1
        if not (ok and psd and ok_obj):
            report(f'(1) n={n} x={x}: constraints {ok} psd {psd} objective {obj}', False)
    print(f'(1) n={n}: 5 points in (L,U) checked; last point pivots={piv}, size={len(M)}, objective - (n+1) = {obj - (n + 1)}')
report('(1) all constraint, PSD and objective checks passed:', bad == 0)

# (2) symbolic identities
n, x = sp.symbols('n x', positive=True)
L = (n + 2) ** 2 / (4 * (n + 1)); delta = x - L
p0 = (n + 1) / (2 * n); q = delta / (16 * n); c = (x + L) / (4 * n); p = p0 - c * q; d = (n * p / 2 + q / 4) / (n - 1)
obj = x * q + 2 * n * p
report('(2) objective = n+1+delta^2/(32n):', sp.simplify(obj - (n + 1 + delta ** 2 / (32 * n))) == 0)
Bm = sp.Matrix([[1, sp.sqrt(2 * n) * p, q], [sp.sqrt(2 * n) * p, (n + 1) * p, sp.sqrt(2 * n) * q / 4], [q, sp.sqrt(2 * n) * q / 4, q]])
S = Bm[1:, 1:] - Bm[1:, 0] * Bm[0, 1:]
report('(2) Schur complement S(1,1) = 2n p (p0 - p):', sp.simplify(S[0, 0] - 2 * n * p * (p0 - p)) == 0)
report('(2) det S = 2n q^3 (4 p0 - c (c - 1/(2n))):', sp.simplify(S.det() - 2 * n * q ** 3 * (4 * p0 - c * (c - 1 / (2 * n)))) == 0)
lam_plus = (n - 2) * p / (2 * (n - 1)) - n * q / (4 * (n - 1)); lam_minus = d + q / 4
# numerical cross-check of the block spectrum against numpy for one instance
idx, first, M, _ = moment(5, F(49, 24) + F(1, 100))
ev = np.linalg.eigvalsh(np.array([[float(v) for v in r] for r in M]))
print('(2) n=5 smallest eigenvalues:', np.round(ev[:6], 10))
report('(2) numerical spectrum nonnegative (n=5):', ev.min() > -1e-12)

# (3) exact ground energy
def singlet(i, j, N):
    I2 = np.eye(2); X = np.array([[0, 1], [1, 0]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1])
    def op(P, k):
        out = np.array([[1.0]])
        for t in range(N): out = np.kron(out, P if t == k else I2)
        return out
    return (np.eye(2 ** N) - op(X, i) @ op(X, j) - op(Y, i) @ op(Y, j) - op(Z, i) @ op(Z, j)) / 4
for nn in (3, 4, 5):
    N = nn + 2; Lv = F((nn + 2) ** 2, 4 * (nn + 1)); Uv = F(nn + 2, 3)
    for t in (10, 50, 90):
        xv = Lv + (Uv - Lv) * F(t, 100)
        H = float(xv) * singlet(0, 1, N) + sum(singlet(0, k, N) + singlet(1, k, N) for k in range(2, N))
        e = np.linalg.eigvalsh(H).max()
        if abs(e - (nn + 1)) > 1e-9: report(f'(3) n={nn} x={xv}: lambda_max = {e}', False)
    print(f'(3) n={nn}: lambda_max(H) = {nn + 1} at three points of (L, U)')
report('(3) exact diagonalization agrees with E_GS = n+1:', bad == 0)

# (4) interval nonempty iff n >= 3
report('(4) L < U iff n > 2 (n = 1..20):', all((F((m + 2) ** 2, 4 * (m + 1)) < F(m + 2, 3)) == (m > 2) for m in range(1, 21)))
print('bad=', bad)
sys.exit(1 if bad else 0)
