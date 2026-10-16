"""Girard–Gour arXiv:1701.02724v3, Eq. (8) (and the fixed-concurrence lower envelope after Eq. (7)).

N(s) = 2 Tr[(s^G)_-],  N2(s) = Tr[(s^G)_-] + 2 Tr[(((s^G)_-)^G)_-],  G = id (x) T,  X_- = (|X| - X)/2.
Eq. (8): f(nu) = (nu/2)(1 + c_nu / sqrt((1 - c_nu)^2 + c_nu^2)) <= N2(s), c_nu = sqrt(2 nu (nu+1)) - nu, nu = N(s).
After (7): N2(s) >= (nu_c/2)(1 + c/sqrt((1-c)^2+c^2)), c = C(s) Wootters, nu_c = sqrt((1-c)^2+c^2) - (1-c).
Also checked: Lemma 2 (lambda <= r^2/(1+2r-r^2), r = Schmidt ratio of the negative eigenvector), Lemma 3 formula
N2 = nu (1 + 2r/(1+r^2))/2, equality family sigma_c, and the Lemma 5 attainment construction.
"""
import math, sys
import numpy as np
from scipy.optimize import minimize

rng = np.random.default_rng(20261017)


def PT(X):
    return X.reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)


def negpart(X):
    w, V = np.linalg.eigh((X + X.conj().T) / 2)
    wn = np.where(w < 0, -w, 0.0)
    return V @ np.diag(wn) @ V.conj().T


def N(s): return 2 * np.trace(negpart(PT(s))).real
def N2(s):
    m = negpart(PT(s)); return np.trace(m).real + 2 * np.trace(negpart(PT(m))).real


def f(nu):
    if nu <= 0: return 0.0
    c = math.sqrt(2 * nu * (nu + 1)) - nu
    return nu / 2 * (1 + c / math.sqrt((1 - c) ** 2 + c ** 2))


def conc(s):
    sy = np.array([[0, -1j], [1j, 0]]); Y = np.kron(sy, sy)
    R = s @ Y @ s.conj() @ Y
    ev = np.sqrt(np.clip(np.sort(np.linalg.eigvals(R).real)[::-1], 0, None))
    return max(0.0, ev[0] - ev[1] - ev[2] - ev[3])


def g_conc(c):
    if c <= 0: return 0.0
    q = math.sqrt((1 - c) ** 2 + c ** 2); nu = q - (1 - c)
    return nu / 2 * (1 + c / q)


def rand_state(rank):
    A = rng.normal(size=(4, rank)) + 1j * rng.normal(size=(4, rank)); s = A @ A.conj().T
    return s / np.trace(s).real


def param_state(x):
    A = (x[:16] + 1j * x[16:]).reshape(4, 4); s = A @ A.conj().T
    return s / np.trace(s).real


def schmidt_ratio(v):
    M = v.reshape(2, 2); sv = np.linalg.svd(M, compute_uv=False); return sv[1] / sv[0]


def main():
    fails = []; worst8 = math.inf; worstc = math.inf; worstL = math.inf; lem3 = 0.0; count = 0
    for k in range(40000):
        s = rand_state(int(rng.integers(1, 5)))
        nu = N(s)
        if nu < 1e-9: continue
        count += 1
        n2 = N2(s)
        worst8 = min(worst8, n2 - f(nu)); worstc = min(worstc, n2 - g_conc(conc(s)))
        w, V = np.linalg.eigh(PT(s)); lam = -w[0]; r = schmidt_ratio(V[:, 0])
        worstL = min(worstL, r ** 2 / (1 + 2 * r - r ** 2) - lam)
        lem3 = max(lem3, abs(n2 - nu * (1 + 2 * r / (1 + r ** 2)) / 2))
    print(f"random states with N > 0: {count}; min N2 - f(N) = {worst8:.3e}; min N2 - g(C) = {worstc:.3e}")
    print(f"Lemma 2 slack min = {worstL:.3e}; Lemma 3 formula max error = {lem3:.2e}")
    if worst8 < -1e-9 or worstc < -1e-9 or worstL < -1e-9 or lem3 > 1e-9: fails.append("random")
    # adversarial: minimize N2 - f(N) over all states, many restarts
    best = math.inf
    for _ in range(60):
        x0 = rng.normal(size=32)
        obj = lambda x: (lambda s: (N2(s) - f(N(s))) if N(s) > 1e-6 else 1.0)(param_state(x))
        res = minimize(obj, x0, method="Nelder-Mead", options={"maxiter": 6000, "xatol": 1e-10, "fatol": 1e-13})
        best = min(best, res.fun)
    print(f"adversarial min of N2 - f(N) over 60 restarts: {best:.3e} (expected >= 0, approaching 0)")
    if best < -1e-7: fails.append(("adversarial", best))
    # equality family sigma_c and Lemma 5 attainment
    phi = np.array([1, 0, 0, 1]) / math.sqrt(2); e10 = np.zeros(4); e10[2] = 1
    eqerr = 0.0
    for c in np.linspace(0.05, 1, 20):
        s = c * np.outer(phi, phi) + (1 - c) * np.outer(e10, e10)
        eqerr = max(eqerr, abs(N2(s) - f(N(s))))
    print(f"sigma_c family: max |N2 - f(N)| = {eqerr:.2e}")
    if eqerr > 1e-9: fails.append("equality")
    att = 0.0
    for _ in range(500):
        r = rng.uniform(0.05, 1); lam = rng.uniform(0, 1) * r ** 2 / (1 + 2 * r - r ** 2)
        h = (1 + 2 * lam) / (r + 1 / r + 2); a = h * r - lam; d = h / r - lam
        s = np.array([[a, 0, 0, 0], [0, h, -h, 0], [0, -h, h, 0], [0, 0, 0, d]], float)
        if np.linalg.eigvalsh(s).min() < -1e-12 or abs(np.trace(s) - 1) > 1e-12: fails.append("lemma5-state"); break
        nu = 2 * lam
        att = max(att, abs(N(s) - nu) + abs(N2(s) - nu * (1 + 2 * r / (1 + r ** 2)) / 2))
    print(f"Lemma 5 construction: max error in (N, N2) = {att:.2e}")
    if att > 1e-9: fails.append("lemma5")
    # control: a wrong (too strong) bound must fail: N2 >= N
    viol = sum(1 for _ in range(2000) if (lambda s: N(s) > 1e-6 and N2(s) < N(s) - 1e-9)(rand_state(3)))
    print(f"control N2 >= N fails on {viol}/2000 random rank-3 states (expected many)")
    if viol == 0: fails.append("control")
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
