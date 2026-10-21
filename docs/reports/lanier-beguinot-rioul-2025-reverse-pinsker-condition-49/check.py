"""Check of Eq. (50) of arXiv:2501.14340v1 for the maximal chi-square divergence on qubit pairs.

sigma = diag(a, 1-a), rho = sigma + t X with X the Pauli X matrix and q = a(1-a).  Then
sigma^{-1/2} rho sigma^{-1/2} = I + (t/sqrt q) X, D^max_f = t^2/q for f(x) = (x-1)^2, and the right side of
Eq. (50) is 2 t^2 / sqrt q.  Part (a) is exact rational arithmetic at a = 9/25, t = 12/125.  Part (b) is floating
point: a sweep over (a, t) and a Hilbert-Schmidt sample of 4x4 pairs.  Part (c) samples 4x4 pairs with a Haar
eigenbasis and uniform eigenvalues and reports the share satisfying condition (49), which the source puts above 80%.
"""
from fractions import Fraction as F
import numpy as np

ok = True
# (a) exact instance
a, t = F(9, 25), F(12, 125)
q = a * (1 - a)
sq = F(12, 25)
assert sq * sq == q
c = t / sq
m, M = 1 - c, 1 + c
Dmax = c * c  # Tr sigma^{1/2} (cX)^2 sigma^{1/2} = c^2 Tr sigma
tracenorm = 2 * t
rhs = tracenorm / 2 * ((m - 1) ** 2 / (1 - m) + (M - 1) ** 2 / (M - 1))
detrho = a * (1 - a) - t * t
# (49): rho + sigma - |rho - sigma| = rho + sigma - t I
P = [[2 * a - t, t], [t, 2 * (1 - a) - t]]
detP = P[0][0] * P[1][1] - P[0][1] * P[1][0]
print(f"(a) a={a} t={t}: sqrt q={sq}, m={m}, M={M}, D^max={Dmax}, ||rho-sigma||_1={tracenorm}, rhs(50)={rhs}")
print(f"    det rho={detrho} > 0: {detrho > 0}; (49) matrix diag {P[0][0]}, {P[1][1]}, det {detP}: PSD {P[0][0] > 0 and detP >= 0}")
print(f"    D^max > rhs(50): {Dmax > rhs} (difference {Dmax - rhs})")
ok &= detrho > 0 and Dmax > rhs and P[0][0] > 0 and detP >= 0


def fun(A, f):
    w, V = np.linalg.eigh(A)
    return V @ np.diag([f(v) for v in w]) @ V.conj().T


def evaluate(rho, sigma, f):
    sh = fun(sigma, np.sqrt)
    shi = np.linalg.inv(sh)
    R = shi @ rho @ shi
    w = np.linalg.eigvalsh((R + R.conj().T) / 2)
    m, M = w.min(), w.max()
    Dmax = np.trace(sh @ fun((R + R.conj().T) / 2, f) @ sh).real
    tn = np.abs(np.linalg.eigvalsh(rho - sigma)).sum()
    rhs = tn / 2 * (f(m) / (1 - m) + f(M) / (M - 1))
    c49 = np.linalg.eigvalsh(rho + sigma - fun(rho - sigma, abs)).min()
    return Dmax, rhs, c49


f = lambda v: (v - 1) ** 2
print("(b) qubit sweep, sigma = diag(a,1-a), rho = sigma + tX:")
for aa in [0.6, 0.75, 0.9]:
    for tt in [0.02, 0.1, 0.2]:
        qq = aa * (1 - aa)
        if tt * tt >= qq:
            continue
        s = np.diag([aa, 1 - aa]).astype(complex)
        r = s + tt * np.array([[0, 1], [1, 0]], dtype=complex)
        D, rh, c49 = evaluate(r, s, f)
        print(f"    a={aa} t={tt}: D^max={D:.6f} rhs={rh:.6f} ratio={D / rh:.4f} (49) min eigenvalue {c49:.4f}")
        ok &= D > rh
rng = np.random.default_rng(20250124)


def rand_state(d):
    G = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    W = G @ G.conj().T
    return W / np.trace(W).real


N = 2000
sat = viol = viol_sat = 0
for _ in range(N):
    rho, sigma = rand_state(4), rand_state(4)
    D, rh, c49 = evaluate(rho, sigma, f)
    s49 = c49 >= -1e-10
    sat += s49
    viol += D > rh * (1 + 1e-9)
    viol_sat += s49 and D > rh * (1 + 1e-9)
print(f"(b) {N} Hilbert-Schmidt random 4x4 pairs (seed 20250124): (49) holds {sat}/{N}; (50) fails for D^max_(x-1)^2 {viol}/{N}; both {viol_sat}/{N}")
def rand_state_spread(d):
    U = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))[0]
    w = rng.random(d)
    return U @ np.diag(w / w.sum()) @ U.conj().T


sat = viol_sat = 0
for _ in range(N):
    rho, sigma = rand_state_spread(4), rand_state_spread(4)
    D, rh, c49 = evaluate(rho, sigma, f)
    s49 = c49 >= -1e-10
    sat += s49
    viol_sat += s49 and D > rh * (1 + 1e-9)
print(f"(c) {N} random 4x4 pairs U diag(w/sum w) U*, U Haar, w uniform on [0,1]^4: (49) holds {sat}/{N} "
      f"(the source reports more than 80% for its sampler); (49) holds and (50) fails for D^max_(x-1)^2 {viol_sat}/{N}")

print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
