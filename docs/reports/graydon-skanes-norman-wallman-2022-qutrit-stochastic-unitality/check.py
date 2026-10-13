#!/usr/bin/env python3
"""A stochastic qutrit channel that is not unital.

Source: Graydon, Skanes-Norman, Wallman, "Designing Stochastic Channels", arXiv:2201.07156v1, Section II:
"We leave open the question of whether or not all stochastic channels are unital for d = 3."
A stochastic channel is a CPTP map with J(Phi) col(1_d) = lambda col(1_d) for some 0 < lambda <= 1, where
J(Phi) = (1/d) (Phi (x) id)(col(1_d) col(1_d)^dagger)  (Eqs. (2), (11)).

Family (every d >= 3, 0 < lam < 1): m = d - 1, D = diag(1, -1/m, ..., -1/m), c^2 = 1 - 1/m^2,
Phi = lam id + (1-lam) Psi with Psi(X) = D X D + c^2 sum_{j>=1} E_{0j} X E_{j0}.
Exact rational checks for d = 3..8 and lam in {1/2, 1/3, 9/10}:
  * sum K^dagger K = I (trace preservation; Kraus form gives complete positivity);
  * J(Phi) col(1) = lam col(1) (all error Kraus operators are traceless and Hilbert-Schmidt orthogonal to I);
  * Phi(I) = diag(1 + (1-lam)(m - 1/m), lam + (1-lam)/m^2, ...) != I.
At d = 3, lam = 1/2: Phi(I) = diag(7/4, 5/8, 5/8).
Prints one line per case and ALL_OK.
"""
from fractions import Fraction as Fr
import sys
def fail(msg): print("FAIL", msg); sys.exit(1)
for d in range(3, 9):
    m = d - 1
    for lam in (Fr(1, 2), Fr(1, 3), Fr(9, 10)):
        # Kraus operators as (weight, matrix) with channel X -> sum w K X K^T (all real)
        I = [[Fr(int(i == j)) for j in range(d)] for i in range(d)]
        Dm = [[(Fr(1) if i == 0 else Fr(-1, m)) if i == j else Fr(0) for j in range(d)] for i in range(d)]
        c2 = 1 - Fr(1, m * m)
        E = lambda i, j: [[Fr(int(r == i and s == j)) for s in range(d)] for r in range(d)]
        kraus = [(lam, I), (1 - lam, Dm)] + [((1 - lam) * c2, E(0, j)) for j in range(1, d)]
        def mul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(d)) for j in range(d)] for i in range(d)]
        def T(A): return [list(r) for r in zip(*A)]
        def add(A, B): return [[A[i][j] + B[i][j] for j in range(d)] for i in range(d)]
        Z = [[Fr(0)] * d for _ in range(d)]
        tp = Z; un = Z
        for w, K in kraus:
            tp = add(tp, [[w * x for x in row] for row in mul(T(K), K)])
            un = add(un, [[w * x for x in row] for row in mul(K, T(K))])
        if tp != I: fail(f"trace preservation d={d} lam={lam}")
        # J col(1) = (1/d) sum_w w col(K) tr(K^T): only K = I has nonzero trace
        traces = [sum(K[i][i] for i in range(d)) for _, K in kraus]
        if any(t != 0 for (w, K), t in zip(kraus[1:], traces[1:])): fail("traceless errors")
        lam_eff = kraus[0][0] * traces[0] / d
        if lam_eff != lam: fail(f"eigenvalue d={d} lam={lam}")
        # Hilbert-Schmidt orthogonality of the error operators to I and to each other
        for a in range(1, len(kraus)):
            for b in range(a + 1, len(kraus)):
                if sum(kraus[a][1][i][j] * kraus[b][1][i][j] for i in range(d) for j in range(d)) != 0: fail("orthogonality")
        if un == I: fail(f"unital d={d} lam={lam}")
        print(f"d={d} lam={lam}: stochastic (J col(1) = {lam_eff} col(1)), Phi(I) diagonal = {[str(un[i][i]) for i in range(d)]}")
print("ALL_OK")
