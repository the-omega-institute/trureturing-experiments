"""Poncini-Rasmussen, "Integrability of planar-algebraic models", arXiv:2206.14462v2, Appendix C, Conjecture C.12 (C.66).

Checks (exact, sympy):
 1. Scalar form: for n = 1..40 and k = 0..n-1, the right side of (C.66) with h_n replaced by its eigenvalue n-1-2k equals
    the eigenvalue 2(1+u^2)(1-u)^{2k}(1+u)^{2(n-1-k)} of (C.37).
 2. Literal form: for n = 1..4, the transfer matrix built from the definitions (C.9)-(C.15) (eight-vertex R, K = 1,
    pre-trace operator (C.12), partial trace (C.13)) equals the right side of (C.66) as a matrix, with h_n from (C.40),
    m_j from (C.44) and the double factorial binomial (C.65) with (-1)!! = 1.
"""
import itertools, sympy as sp
u = sp.symbols('u')
def dfac(n):
    if n in (-1, 0): return 1
    r = 1
    while n > 0: r *= n; n -= 2
    return r
def dfbin(a, b): return sp.Rational(dfac(a), dfac(b) * dfac(a - b))
def m(j, x):
    r = 1
    for k in range(j): r *= (x - (j - 1 - 2 * k))
    return r
def rhs_scalar(n, x):
    return 2 * sum((1 + u**2)**(n - k) * (2*u)**k * sum((-1)**j * dfbin(n - k - 2 + 2*j, 2*j) * m(k - 2*j, x) / sp.factorial(k - 2*j)
                   for j in range(k // 2 + 1)) for k in range(n))
# scalar identity: eigenvalue lambda_k = 2(1+u^2)(1-u)^{2k}(1+u)^{2(n-1-k)} at h-eigenvalue mu_k = n-1-2k
ok = True
for n in range(1, 41):
    for k in range(n):
        lam = 2 * (1 + u**2) * (1 - u)**(2*k) * (1 + u)**(2*(n - 1 - k))
        if sp.expand(rhs_scalar(n, n - 1 - 2*k) - lam) != 0:
            ok = False; print('scalar FAIL', n, k)
print('scalar identity (eigenvalues of C.37 at h-eigenvalues n-1-2k) for n=1..40:', ok)
# literal matrix definition for n = 1..4
import numpy as np
def kron(*ms):
    r = sp.Matrix([[1]])
    for M in ms: r = sp.kronecker_product(r, M)
    return r
def E(k, j):
    M = sp.zeros(2, 2); M[k, j] = 1; return M
def Rcoef(j, k, l, mm):  # R^{km}_{jl}(u), indices 0/1
    if j == mm and l == k: return 1
    if j == l and k == mm and j != k: return u
    if j == k and l == mm and j != l: return u
    return 0
def R_i_last(i, N):  # R_{i,N} on N factors (1-based i < N): sum R^{km}_{jl} 1 x e_kj (at i) x 1 x e_ml (at N)
    res = sp.zeros(2**N, 2**N)
    for j, k, l, mm in itertools.product(range(2), repeat=4):
        c = Rcoef(j, k, l, mm)
        if c == 0: continue
        fs = [sp.eye(2)] * N; fs = list(fs); fs[i-1] = E(k, j); fs[N-1] = E(mm, l)
        res += c * kron(*fs)
    return res
def R_last_i(i, N):  # R_{N,i}: sum R^{km}_{jl} 1 x e_ml (at i) x 1 x e_kj (at N)
    res = sp.zeros(2**N, 2**N)
    for j, k, l, mm in itertools.product(range(2), repeat=4):
        c = Rcoef(j, k, l, mm)
        if c == 0: continue
        fs = [sp.eye(2)] * N; fs = list(fs); fs[i-1] = E(mm, l); fs[N-1] = E(k, j)
        res += c * kron(*fs)
    return res
def ptrace_last(M, N):
    d = 2**(N-1); out = sp.zeros(d, d)
    for a in range(d):
        for b in range(d):
            out[a, b] = M[2*a, 2*b] + M[2*a+1, 2*b+1]
    return out
X = sp.Matrix([[0, 1], [1, 0]])
mat_ok = True
for n in range(1, 5):
    N = n + 1
    L = sp.eye(2**N)
    for i in range(n, 0, -1): L = L * R_i_last(i, N)
    for i in range(1, n + 1): L = L * R_last_i(i, N)
    T = sp.expand(ptrace_last(L, N))
    h = sp.zeros(2**n, 2**n)
    for i in range(1, n):
        fs = [sp.eye(2)] * n; fs = list(fs); fs[i-1] = X; fs[i] = X; h += kron(*fs)
    def mpoly(j):
        r = sp.eye(2**n)
        for k in range(j): r = r * (h - (j - 1 - 2*k) * sp.eye(2**n))
        return r
    Z0 = sp.zeros(2**n, 2**n)
    R = 2 * sum((((1 + u**2)**(n - k) * (2*u)**k) * sum((((-1)**j * dfbin(n - k - 2 + 2*j, 2*j) / sp.factorial(k - 2*j)) * mpoly(k - 2*j)
                for j in range(k // 2 + 1)), Z0) for k in range(n)), Z0)
    diff = sp.expand(T - R)
    good = diff == sp.zeros(2**n, 2**n)
    mat_ok &= good
    print('n', n, 'literal transfer matrix equals C.66:', good)
print('RESULT', 'PASS' if (ok and mat_ok) else 'FAIL')
raise SystemExit(0 if (ok and mat_ok) else 1)
