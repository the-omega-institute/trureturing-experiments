#!/usr/bin/env python3
"""Exact checks behind the proof that no symmetry of the set of quantum channels reverses time.

Source: Chiribella, Aurell, Zyczkowski, "Symmetries of quantum evolutions", arXiv:2101.04962v5;
Phys. Rev. Research 3, 033028 (2021), Section V.B, the conjecture after Theorem 4.

Part I (d >= 3): for P with columns e1, e1, e3, ..., ed, r = P 1 and Q = P - (r - 1) 1^T / d
must equal G(P)^T, a column-stochastic matrix; Q has row and column sums 1 and Q[0][2] = -1/d < 0.
Part II (d = 2): the amplitude-damping channel with Kraus A0 = [[1,0],[0,1/2]], A1 = [[0,sqrt3/2],[0,0]]
is trace preserving with Bloch data K = diag(1/2,1/2,1/4), t = (0,0,3/4); the bound |w_j| <= 3/16
forces |w|^2 <= 27/256 < 9/16 = |t|^2.
Part III: tau o Ad_U o tau = Ad_{conj U}, so Ad_{U^T} maps to Ad_{U^dagger}; checked on random unitaries.
Prints one line per check and ALL_OK.
"""
from fractions import Fraction as Fr
import cmath, math, random, sys

def fail(m): print("FAIL", m); sys.exit(1)

# Part I
for d in range(3, 9):
    P = [[Fr(0)] * d for _ in range(d)]
    cols = [0, 0] + list(range(2, d))
    for j, i in enumerate(cols): P[i][j] = Fr(1)
    r = [sum(P[i]) for i in range(d)]
    Q = [[P[i][j] - (r[i] - 1) / d for j in range(d)] for i in range(d)]
    if any(sum(Q[i]) != 1 for i in range(d)): fail(f"row sums d={d}")
    if any(sum(Q[i][j] for i in range(d)) != 1 for j in range(d)): fail(f"col sums d={d}")
    if Q[0][2] != Fr(-1, d): fail(f"entry d={d}")
    print(f"PART_I d={d}: r={[int(x) for x in r]} Q[1][3]={Q[0][2]} (forced G(P)[3][1] < 0)")

# Part II: amplitude damping, exact arithmetic with s = sqrt(3)/2 handled symbolically via s^2 = 3/4
s2 = Fr(3, 4)
# A0^dag A0 + A1^dag A1 = diag(1, 1/4) + diag(0, s^2) = I
if (Fr(1), Fr(1, 4) + s2) != (Fr(1), Fr(1)): fail("trace preservation")
# Bloch data: rho_v = (I + v.sigma)/2 -> output; t = image of v=0, K diag from output on axes
# channel on diag(a, b): a' = a + s^2 b, b' = b/4 ; on off-diagonal c: c' = c/2
t3 = (Fr(1, 2) + s2 * Fr(1, 2)) - (Fr(1, 2) * Fr(1, 4))   # z-component of output on maximally mixed
K = (Fr(1, 2), Fr(1, 2), Fr(1, 4))
if t3 != Fr(3, 4): fail("t")
print(f"PART_II amplitude damping: K=diag{tuple(str(k) for k in K)} t=(0,0,{t3}); 13/16 + w_j <= 1 gives |w_j| <= 3/16; 3*(3/16)^2 = {3*Fr(3,16)**2} < {t3**2}")
if not 3 * Fr(3, 16) ** 2 < t3 ** 2: fail("contradiction")
if Fr(1, 2) ** 2 + t3 ** 2 != Fr(13, 16): fail("13/16")

# Part III: transpose conjugation identity on random 3x3 unitaries (floating point, tolerance 1e-12)
def mat(n, f): return [[f(i, j) for j in range(n)] for i in range(n)]
def mul(A, B): n = len(A); return mat(n, lambda i, j: sum(A[i][k] * B[k][j] for k in range(n)))
def T(A): n = len(A); return mat(n, lambda i, j: A[j][i])
def H(A): n = len(A); return mat(n, lambda i, j: A[j][i].conjugate())
def conjm(A): n = len(A); return mat(n, lambda i, j: A[i][j].conjugate())
random.seed(1)
for trial in range(5):
    n = 3
    M = mat(n, lambda i, j: complex(random.gauss(0, 1), random.gauss(0, 1)))
    # Gram-Schmidt to a unitary
    cols = [[M[i][j] for i in range(n)] for j in range(n)]; U = []
    for c in cols:
        for u in U:
            ip = sum(u[i].conjugate() * c[i] for i in range(n)); c = [c[i] - ip * u[i] for i in range(n)]
        nr = math.sqrt(sum(abs(x) ** 2 for x in c)); U.append([x / nr for x in c])
    U = T(U)
    X = mat(n, lambda i, j: complex(random.gauss(0, 1), random.gauss(0, 1)))
    lhs = T(mul(mul(U, T(X)), H(U)))          # tau o Ad_U o tau (X)
    rhs = mul(mul(conjm(U), X), T(U))         # Ad_{conj U}(X)
    if max(abs(lhs[i][j] - rhs[i][j]) for i in range(n) for j in range(n)) > 1e-12: fail("part III")
print("PART_III tau o Ad_U o tau = Ad_{conj U} on 5 random unitaries")
print("ALL_OK")
