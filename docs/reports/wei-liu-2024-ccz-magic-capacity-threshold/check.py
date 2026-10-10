# Check of a refutation of Wei–Liu (arXiv:2410.21215v1, Sec. VII): "We conjecture that the magic capacity threshold
# for CCZ under local depolarizing noise is 1/3."  Witness: a Bell-type stabilizer input with CNOT flags on the reference.
import numpy as np, itertools, random
from fractions import Fraction
random.seed(7); np.random.seed(7)
I2 = np.eye(2); H = np.array([[1,1],[1,-1]])/np.sqrt(2); S = np.diag([1,1j])
X = np.array([[0,1],[1,0]]); Z = np.diag([1,-1]); Y = 1j*X@Z
def kron(*ms):
    r = np.array([[1]])
    for m in ms: r = np.kron(r, m)
    return r
def op1(m, i, n):
    return kron(*[m if k == i else I2 for k in range(n)])
def cnot(c, t, n):
    d = 2**n; M = np.zeros((d, d))
    for x in range(d):
        bits = [(x >> (n-1-k)) & 1 for k in range(n)]
        if bits[c]: bits[t] ^= 1
        y = sum(b << (n-1-k) for k, b in enumerate(bits)); M[y, x] = 1
    return M
CCZ = np.diag([1,1,1,1,1,1,1,-1])
plus3 = np.ones(8)/np.sqrt(8); psi = CCZ @ plus3
# 1. all 3-qubit pure stabilizer states by BFS over the Clifford generators
def canon(v):
    k = np.argmax(np.abs(v) > 1e-9); v = v * (abs(v[k]) / v[k])
    return tuple(np.round(v.real, 9)) + tuple(np.round(v.imag, 9))
gens = [op1(H,i,3) for i in range(3)] + [op1(S,i,3) for i in range(3)] + [cnot(c,t,3) for c in range(3) for t in range(3) if c != t]
start = np.zeros(8, complex); start[0] = 1
seen = {canon(start): start}; frontier = [start]
while frontier:
    nxt = []
    for v in frontier:
        for g in gens:
            w = g @ v; k = canon(w)
            if k not in seen: seen[k] = w; nxt.append(w)
    frontier = nxt
stab3 = list(seen.values())
assert len(stab3) == 1080, len(stab3)
F = max(abs(np.vdot(s, psi))**2 for s in stab3)
assert abs(F - 9/16) < 1e-12
print('3-qubit stabilizer states:', len(stab3), '; max |<s|CCZ+++>|^2 =', F, '(9/16)')
# 2. the noisy channel on A with reference B, and the witness
n = 6
Omega = np.zeros(64, complex)
for xb in range(8): Omega[(xb << 3) | xb] = 1/np.sqrt(8)       # qubits A1A2A3 B1B2B3
U = kron(CCZ, np.eye(8))
def depol_A(rho, lam):
    for i in range(3):
        paulis = [op1(P, i, 6) for P in (X, Y, Z)]
        rho = (1 - 3*lam/4) * rho + (lam/4) * sum(P @ rho @ P.conj().T for P in paulis)
    return rho
V = cnot(0,3,6) @ cnot(1,4,6) @ cnot(2,5,6)
P000 = np.zeros((8,8)); P000[0,0] = 1
W = V.conj().T @ np.kron(9/16*np.eye(8) - np.outer(psi, psi.conj()), P000) @ V
for lam in [0, 0.1, 1/3, 0.4, 0.5, 0.51, 0.52]:
    rho = depol_A(U @ np.outer(Omega, Omega.conj()) @ U.conj().T, lam)
    w = np.trace(W @ rho).real
    formula = (9/16)*(1-lam/2)**3 - (1-3*lam/4)**3
    assert abs(w - formula) < 1e-12, (lam, w, formula)
    print(f'lambda={lam:.4f}: tr(W rho) = {w:.6f}  (formula {formula:.6f})')
rho = depol_A(U @ np.outer(Omega, Omega.conj()) @ U.conj().T, 0.5)
assert abs(np.trace(W @ rho).real - (-7/1024)) < 1e-12
print('lambda = 1/2: tr(W rho) = -7/1024 exactly (float check)')
# exact sign on [0, 1/2] of w(lam) = (9/16)(1-lam/2)^3 - (1-3lam/4)^3 via rationals on a fine grid plus monotone factor
for k in range(0, 501):
    l = Fraction(k, 1000)
    assert Fraction(9,16)*(1-l/2)**3 - (1-3*l/4)**3 < 0
print('w(lambda) < 0 at every lambda = k/1000 in [0, 1/2] (exact rationals)')
# 3. W >= 0 on random 6-qubit stabilizer states (sanity check of the measurement-closure argument)
gens6 = [op1(H,i,6) for i in range(6)] + [op1(S,i,6) for i in range(6)] + [cnot(c,t,6) for c in range(6) for t in range(6) if c != t]
mn = 1e9
for trial in range(3000):
    v = np.zeros(64, complex); v[0] = 1
    for _ in range(80): v = random.choice(gens6) @ v
    mn = min(mn, np.vdot(v, W @ v).real)
assert mn > -1e-12
print('min <phi|W|phi> over 3000 random 6-qubit stabilizer states =', round(mn, 12))
print('ALL_OK')
# 4. the other composition order (noise before the gate): same witness
for lam in [1/3, 0.5]:
    rho2 = U @ depol_A(np.outer(Omega, Omega.conj()), lam) @ U.conj().T
    w2 = np.trace(W @ rho2).real
    print(f'noise before CCZ, lambda={lam:.4f}: tr(W rho) = {w2:.6f}')
    assert w2 < 0
print('ALL_OK_BOTH_ORDERS')
