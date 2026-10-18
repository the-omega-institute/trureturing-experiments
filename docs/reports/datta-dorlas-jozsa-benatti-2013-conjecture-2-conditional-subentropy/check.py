"""Datta-Dorlas-Jozsa-Benatti (arXiv:1310.1312) Conjecture 2: Q(X|B) <= -ln p_guess(X|B).

Checks, with Q from Equation (1) (distinct eigenvalues, 60-digit arithmetic):
 1. Equation (4) agrees with Equation (1) (n = 2, 3, quadrature).
 2. Random c-q states: 0 <= Q(rho_XB) - Q(rho_B) <= 1 - P_up <= -ln P_up, where P_up >= p_guess is a
    certified dual bound tr(sigma) with sigma >= A_x for every x (so the checked inequality implies the
    proposed bound Q(X|B) <= 1 - p_guess and hence the conjecture on these states).
 3. Lemma 1: Q(r) <= 1 - max_j r_j for random probability vectors r, and f'(0+) = 1 for
    f(t) = Q(1 - t, t/(m-1), ..., t/(m-1)).
 4. Lemma 2 by Monte Carlo on the complex sphere: Q(Delta R) - Q(R) = D_KL(q_R || q_{Delta R}).
"""
import numpy as np, mpmath as mp, cvxpy as cp
mp.mp.dps = 60
rng = np.random.default_rng(20261018)
ok = True

def Qspec(lam, tol=1e-13):
    lam = [l if isinstance(l, mp.mpf) else mp.mpf(float(l)) for l in lam if l > tol]
    N = len(lam)
    return -sum(li ** N * mp.log(li) / mp.fprod([li - lj for j, lj in enumerate(lam) if j != i]) for i, li in enumerate(lam))
eta = lambda y: -y * mp.log(y) if y > 0 else mp.mpf(0)
Cn = lambda n: sum(mp.mpf(1) / k for k in range(2, n + 1))
# 1. Equation (4)
mp.mp.dps = 30
for lam in ([mp.mpf('0.7'), mp.mpf('0.3')], [mp.mpf('0.5'), mp.mpf('0.3'), mp.mpf('0.2')], [mp.mpf('0.61'), mp.mpf('0.27'), mp.mpf('0.12')]):
    n = len(lam)
    if n == 2:
        q4 = 2 * mp.quad(lambda t: eta(lam[0] * t + lam[1] * (1 - t)), [0, 1]) - Cn(2)
    else:
        q4 = 6 * mp.quad(lambda s: mp.quad(lambda t: eta(lam[0] * s + lam[1] * t + lam[2] * (1 - s - t)), [0, 1 - s]), [0, 1]) - Cn(3)
    err = abs(q4 - Qspec(lam)); print('eq (4) vs eq (1)', [mp.nstr(l, 3) for l in lam], mp.nstr(err, 3)); ok &= err < 1e-20
mp.mp.dps = 60
# 2. Random c-q states
def rand_psd(d, r):
    G = rng.normal(size=(d, r)) + 1j * rng.normal(size=(d, r)); return G @ G.conj().T
def p_guess_bounds(As):
    """Primal and dual SDPs for p_guess; the dual sigma is repaired to an exact feasible point."""
    d = As[0].shape[0]; m = len(As)
    E = [cp.Variable((d, d), hermitian=True) for _ in range(m)]
    prim = cp.Problem(cp.Maximize(cp.real(sum(cp.trace(E[x] @ As[x]) for x in range(m)))),
                      [e >> 0 for e in E] + [sum(E) == np.eye(d)])
    prim.solve(solver=cp.CLARABEL)
    # a feasible POVM from the solver output: clip negative parts, then renormalize
    Ev = []
    for e in E:
        w, V = np.linalg.eigh((e.value + e.value.conj().T) / 2); Ev.append((V * np.clip(w, 0, None)) @ V.conj().T)
    T = sum(Ev); w, V = np.linalg.eigh(T); Ti = (V / np.sqrt(w)) @ V.conj().T; Ev = [Ti @ e @ Ti for e in Ev]
    p_low = sum(np.trace(e @ A).real for e, A in zip(Ev, As))
    sig = cp.Variable((d, d), hermitian=True)
    dual = cp.Problem(cp.Minimize(cp.real(cp.trace(sig))), [sig - A >> 0 for A in As])
    dual.solve(solver=cp.CLARABEL)
    S = (sig.value + sig.value.conj().T) / 2
    shift = max(0.0, max(np.linalg.eigvalsh(A - S).max() for A in As))
    return p_low, np.trace(S).real + shift * d
worst = -1.0; ntest = 0; max_gap = 0.0
for trial in range(120):
    d = int(rng.integers(2, 5)); m = int(rng.integers(2, 5))
    ranks = [int(rng.integers(1, d + 1)) for _ in range(m)]
    As = [rand_psd(d, r) for r in ranks]
    if trial % 3 == 0:  # nearly distinguishable ensembles stress the bound
        U = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))[0]
        As = [A + 6 * np.outer(U[:, x % d], U[:, x % d].conj()) for x, A in enumerate(As)]
    t = sum(np.trace(A).real for A in As); As = [A / t for A in As]
    lamXB = np.concatenate([np.linalg.eigvalsh(A) for A in As]); lamB = np.linalg.eigvalsh(sum(As))
    gap = Qspec(lamXB) - Qspec(lamB)
    p_low, p_up = p_guess_bounds(As)
    ok &= gap >= -1e-30 and p_low <= p_up + 1e-12; max_gap = max(max_gap, p_up - p_low)
    margin = float(gap - (1 - p_up)); worst = max(worst, margin); ntest += 1
    ok &= margin < 0 and (1 - p_up) <= -np.log(p_up) + 1e-15
print('c-q states tested:', ntest, 'max of Q(X|B) - (1 - P_up):', worst, ' max duality gap:', max_gap)
# 3. Lemma 1
w1 = -1.0
for trial in range(400):
    m = int(rng.integers(2, 8)); r = rng.dirichlet(np.ones(m) * rng.choice([0.3, 1, 3]))
    w1 = max(w1, float(Qspec(r) - (1 - r.max())))
print('Lemma 1: max of Q(r) - (1 - max r) over 400 vectors:', w1); ok &= w1 < 0
mp.mp.dps = 250
for m in (2, 3, 5, 8):
    t = mp.mpf('1e-20'); v = [1 - t] + [t / (m - 1)] * (m - 1)
    # perturb ties so that Equation (1) applies; the perturbation is O(1e-30)
    v = [v[0]] + [v[k] * (1 + mp.mpf(k) * mp.mpf('1e-40')) for k in range(1, m)]
    s = sum(v); v = [vi / s for vi in v]
    der = Qspec(v, tol=0) / t; print('f\'(0+) at m =', m, ':', mp.nstr(der, 10)); ok &= abs(der - 1) < 1e-6
mp.mp.dps = 60
# 4. Lemma 2 by Monte Carlo
N = 4; R = rand_psd(N, N); R /= np.trace(R).real
P = [np.diag([1, 1, 0, 0]), np.diag([0, 0, 1, 0]), np.diag([0, 0, 0, 1])]
DR = sum(Pk @ R @ Pk for Pk in P)
g = rng.normal(size=(400000, N)) + 1j * rng.normal(size=(400000, N)); psi = g / np.linalg.norm(g, axis=1, keepdims=True)
qR = N * np.einsum('ij,jk,ik->i', psi.conj(), R, psi).real; qD = N * np.einsum('ij,jk,ik->i', psi.conj(), DR, psi).real
kl = np.mean(qR * np.log(qR / qD)); lhs = float(Qspec(np.linalg.eigvalsh(DR)) - Qspec(np.linalg.eigvalsh(R)))
se = np.std(qR * np.log(qR / qD)) / np.sqrt(len(qR))
print('Lemma 2: Q(Delta R) - Q(R) =', round(lhs, 5), ' Monte Carlo KL =', round(kl, 5), '+-', round(se, 5)); ok &= abs(lhs - kl) < 5 * se
print('RESULT', 'PASS' if ok else 'FAIL')
