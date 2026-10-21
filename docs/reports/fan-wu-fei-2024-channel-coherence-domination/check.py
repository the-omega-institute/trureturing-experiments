"""Check of the conjecture C_{alpha,1}(phi) >= Ctilde_{alpha,1}(phi) of Fan, Wu and Fei (QIP 23, 100 (2024); arXiv:2506.10796v1).

D_{alpha,1}(rho, sigma) = ((Tr sigma^{(1-alpha)/2} rho^alpha sigma^{(1-alpha)/2})^{1/alpha} - 1)/(alpha - 1) (Eqs. (1)-(2), z = 1).
C_{alpha,1}(phi) = min over incoherent (diagonal) states sigma on A x B of D(M_phi, sigma), M_phi the normalized Choi
matrix (Eq. (7)); Ctilde_{alpha,1}(phi) = sup_rho D(phi(Delta rho), Delta(phi rho)) (Eq. (10)).
(a) Convention check: for the phase-flip channel, numerical minimization over diagonal states reproduces the source's
Eq. (13).  (b) The qubit channel phi(rho) = rho_00 sigma_r + rho_11 sigma_s, sigma_t = (I + t X)/2, 0 < s < r < 1:
exact values at alpha = 2, r = 3/4, s = 1/2, and a numerical sweep over alpha in (0,1) u (1,2].
"""
from fractions import Fraction as Fr
import math
import numpy as np

rng = np.random.default_rng(11)


def mpow(A, p):
    w, V = np.linalg.eigh((A + A.conj().T) / 2)
    w = np.clip(w, 0, None)
    return V @ np.diag([x ** p if x > 1e-300 else 0.0 for x in w]) @ V.conj().T


def D(alpha, rho, sigma):
    s = mpow(sigma, (1 - alpha) / 2)
    f = np.trace(s @ mpow(rho, alpha) @ s).real
    return (f ** (1 / alpha) - 1) / (alpha - 1)


def choi(phi, dA, dB):
    M = np.zeros((dA * dB, dA * dB), complex)
    for i in range(dA):
        for j in range(dA):
            E = np.zeros((dA, dA), complex); E[i, j] = 1
            M[i * dB:(i + 1) * dB, j * dB:(j + 1) * dB] = phi(E) / dA
    return M


def C_min(alpha, M, trials=2000):
    """Value at the minimizer sigma_i ∝ <i|M^alpha|i>^(1/alpha) (Corollary 2 of the source's [6]), and the least value
    over random diagonal states (which must not be smaller)."""
    d = np.diag(mpow(M, alpha)).real
    w = d ** (1 / alpha); w = w / w.sum()
    at_opt = D(alpha, M, np.diag(w).astype(complex))
    rnd = min(D(alpha, M, np.diag(rng.dirichlet(np.ones(len(d)))).astype(complex)) for _ in range(trials))
    return at_opt, rnd


def eq8(alpha, M):
    d = np.diag(mpow(M, alpha)).real
    return (np.sum(d ** (1 / alpha)) - 1) / (alpha - 1)


ok = True
X = np.array([[0, 1], [1, 0]], complex); Z = np.array([[1, 0], [0, -1]], complex)
for alpha in (0.5, 2.0):
    for p in (0.2, 0.7):
        phi = lambda r, p=p: p * r + (1 - p) * Z @ r @ Z
        M = choi(phi, 2, 2)
        src = (2 ** (1 - 1 / alpha) * (p ** alpha + (1 - p) ** alpha) ** (1 / alpha) - 1) / (alpha - 1)
        (num, rnd), e8 = C_min(alpha, M), eq8(alpha, M)
        print(f"(a) phase flip p={p} alpha={alpha}: source Eq. (13) = {src:.6f}, Eq. (8) = {e8:.6f}, D at the minimizing diagonal state = {num:.6f}, least over random diagonal states = {rnd:.6f}")
        ok &= abs(src - e8) < 1e-9 and abs(num - e8) < 1e-9 and rnd >= e8 - 1e-9
sig = lambda t: (np.eye(2) + t * X) / 2
r, s = 0.75, 0.5
phi = lambda rho: rho[0, 0] * sig(r) + rho[1, 1] * sig(s)
M = choi(phi, 2, 2)
# exact alpha = 2: C <= D(M, I/4) = sqrt(4 Tr M^2) - 1 = sqrt(45/32) - 1; Ctilde >= D(sigma_r, I/2) = sqrt(1 + r^2) - 1 = 1/4
trM2 = Fr(1, 4) ** 2 * 4 + 2 * (Fr(3, 16) ** 2 + Fr(1, 8) ** 2)
print(f"(b) exact: Tr M^2 = {trM2} (= 45/128: {trM2 == Fr(45, 128)}); C_2,1 <= sqrt(45/32) - 1 = {math.sqrt(45 / 32) - 1:.6f}; "
      f"Ctilde_2,1 >= sqrt(1 + (3/4)^2) - 1 = {math.sqrt(1 + 9 / 16) - 1:.6f}; 45/32 < 25/16: {Fr(45, 32) < Fr(25, 16)}")
ok &= trM2 == Fr(45, 128) and Fr(45, 32) < Fr(25, 16)
print(f"    numerically: D(M, I/4) = {D(2, M, np.eye(4) / 4):.6f}, Eq. (8) value of C_2,1 = {eq8(2, M):.6f}, D(sigma_r, I/2) = {D(2, sig(r), np.eye(2) / 2):.6f}")
for alpha in (0.1, 0.3, 0.5, 0.8, 0.95, 1.05, 1.3, 1.6, 2.0):
    Cup = D(alpha, M, np.eye(4) / 4)
    Clow = D(alpha, sig(r), np.eye(2) / 2)
    e8 = eq8(alpha, M)
    print(f"    alpha={alpha}: C <= D(M, I/4) = {Cup:.6f} (Eq. (8): {e8:.6f}) < Ctilde >= {Clow:.6f}: {Cup < Clow}")
    ok &= Cup < Clow and e8 <= Cup + 1e-12
print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
