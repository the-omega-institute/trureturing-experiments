"""Petrov-Scott, Forum Math. Sigma 14 (2026) e15 (arXiv:2412.21126v2), Conjecture 6.3: the alternative
q-Charlier specialization is Fibonacci positive for some rho, q. Checks the proof for rho = q = 1/4."""
import sys, itertools
from fractions import Fraction as Fr
import sympy as sp

ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and bool(cond)

r, q = sp.symbols('rho q', positive=True)
def xs(k): return q**(k-1)*(1 + r*q**(k-2) + r*q**(k-1) - r*q**(2*k-2)) / ((1 + r*q**(2*k-3))*(1 + r*q**(2*k-1)))
def ys(k): return r*q**(3*k-2)*(1-q**k)*(1 + r*q**(k-1)) / ((1 + r*q**(2*k))*(1 + r*q**(2*k-1))**2*(1 + r*q**(2*k-2)))
def cs(k): return q**(k-1)*(1 + r*q**(k-1)) / ((1 + r*q**(2*k-1))*(1 + r*q**(2*k-2)))
def ts(k): return 0 if k == 0 else r*q**(k-1)*(1-q**k)*(1 + r*q**(2*k+1)) / ((1 + r*q**k)*(1 + r*q**(2*k-1)))
norm = all(sp.simplify(xs(k) - cs(k)*(1 + ts(k-1))) == 0 and sp.simplify(ys(k) - cs(k)*cs(k+1)*ts(k)) == 0 for k in range(1, 9))
check("Lemma 1 (symbolic, k = 1..8): x_k = c_k (1 + t_{k-1}), y_k = c_k c_{k+1} t_k, t_0 = 0", norm)

R = Q = Fr(1, 4)
def ev(e):
    v = sp.Rational(sp.sympify(e).subs({r: sp.Rational(1, 4), q: sp.Rational(1, 4)}))
    return Fr(int(v.p), int(v.q))
K = 60
X = [None] + [ev(xs(k)) for k in range(1, K)]
Y = [None] + [ev(ys(k)) for k in range(1, K)]
T = [ev(ts(k)) if k else Fr(0) for k in range(K)]

def det(M):                       # exact Bareiss-free elimination over Q
    M = [row[:] for row in M]; n = len(M); d = Fr(1)
    for i in range(n):
        p = next((j for j in range(i, n) if M[j][i] != 0), None)
        if p is None: return Fr(0)
        if p != i: M[i], M[p] = M[p], M[i]; d = -d
        d *= M[i][i]
        for j in range(i + 1, n):
            f = M[j][i] / M[i][i]
            for c in range(i, n): M[j][c] -= f * M[i][c]
    return d

def A(l, m):                      # A_l(x+m | y+m)
    return det([[X[i+1+m] if i == j else Y[i+1+m] if j == i+1 else Fr(1) if i == j+1 else Fr(0) for j in range(l)] for i in range(l)])
def B(k, m):                      # B_k(x+m | y+m), (k+1) x (k+1)
    n = k + 1
    def e(i, j):
        if i == 0 and j == 0: return Y[1+m]
        if i == 0 and j == 1: return X[1+m]*Y[2+m]
        if i == j: return X[i+2+m]
        if j == i + 1: return Y[i+2+m]
        if i == j + 1: return Fr(1)
        return Fr(0)
    return det([[e(i, j) for j in range(n)] for i in range(n)])

def s(w):                         # Definition 2.5
    k = 0
    for idx, a in enumerate(w):
        if a == 1: k += 1
        else:
            u = w[idx+1:]; return B(k, sum(u)) * s(u)
    return A(k, 0)

N = 12
words = [w for n in range(0, N + 1) for w in itertools.product((1, 2), repeat=n) if sum(w) <= N]
vals = {w: s(list(w)) for w in words}
check(f"s_w > 0 for all {len(words)} words with |w| <= {N} (exact determinants)", all(v > 0 for v in vals.values()))
check("all t_k in (0, 1/4] and t_{k+1} <= (125/256) t_k for k = 1..58", all(0 < T[k] <= Fr(1, 4) for k in range(1, K)) and all(T[k+1] <= Fr(125, 256)*T[k] for k in range(1, K-1)))
delta = Fr(143, 768)
# normalized B: B_k(m) = (prod c) * E_k(m); check E_k(m) = B_k(m)/prod_{j=1}^{k+2} c_{m+j}... use the continuant form directly
def D(n, m): return 1 + sum(eval_prod(m, rr) for rr in range(1, n + 1))
def eval_prod(m, rr):
    p = Fr(1)
    for j in range(rr): p *= T[m + j]
    return p
def E(k, m):
    if k == 0: return T[m+1]
    return T[m+1] - (1 + T[m] - T[m+1]) * T[m+2] * D(k - 1, m + 3)
check("Lemma 3 bound E_k(m) >= (143/768) t_{m+1} for k, m <= 20", all(E(k, m) >= delta * T[m+1] for k in range(21) for m in range(21)))
C = [None] + [ev(cs(k)) for k in range(1, K)]
def prodc(lo, hi):
    p = Fr(1)
    for j in range(lo, hi + 1): p *= C[j]
    return p
consistent = all(B(k, m) == E(k, m) * prodc(m + 1, m + k + 2) for k in range(0, 8) for m in range(0, 8)) and all(A(l, m) == D(l, m) * prodc(m + 1, m + l) for l in range(0, 8) for m in range(0, 8))
check("Lemma 2 (k, l, m < 8): B_k(m) = E_k(m) prod_{j=m+1}^{m+k+2} c_j and A_l(m) = D_l(m) prod_{j=m+1}^{m+l} c_j", consistent)
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
