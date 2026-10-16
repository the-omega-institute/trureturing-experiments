"""Zhang–Wu–Fei arXiv:1210.4720v4 (PLA 377, 2013), Section 4, conjecture (4.2): counterexample families.

G(rho, sigma; U, V, W) = D(U rho U^dag || sigma) - D(V rho_A V^dag || sigma_A) - D(W rho_B W^dag || sigma_B).
By the source's (4.1) (Zhang–Fei), the supremum over U, V, W of G equals
  H(lambda_down(rho) || lambda_up(sigma)) - H(lambda_down(rho_A) || lambda_down(sigma_A)) - H(lambda_down(rho_B) || lambda_down(sigma_B)),
which is the source's (4.3). Conjecture (4.2) says this supremum is >= 0 for every full-rank sigma.

1. All-dimensional family: rho = I/d, sigma_t = t I + (1 - d t)|00><00|, d = ab, 0 < t < 1/d.
   Direct matrix evaluation (log via eigendecomposition) vs the closed form Delta_ab(t); Delta_ab(t) -> -inf.
2. Two-qubit nondegenerate family: rho_eps = I/4 + eps H (H with entries 3 at 00-11, diag(0,1,-1,0)),
   sigma_eta = diag(61/64, 1/64 - eta, 1/64, 1/64 + eta): sup over unitaries < -1/256 on a grid of
   0 < eps <= 1/1024, 0 < eta <= 1/4096; [rho, sigma] != 0; all spectra simple.
3. Random unitaries never beat the spectral supremum (sanity of the (4.1) reduction), and a control:
   for product sigma = sigma_A (x) sigma_B and rho = I/d the supremum is >= 0 (equality), as additivity predicts.
"""
import math, sys
import numpy as np

rng = np.random.default_rng(20261017)


def logm_h(X):
    w, V = np.linalg.eigh((X + X.conj().T) / 2)
    return V @ np.diag(np.log(w)) @ V.conj().T


def D(r, s):
    w = np.linalg.eigvalsh((r + r.conj().T) / 2)
    lr = logm_h(r) if w.min() > 0 else None
    if lr is None:
        raise ValueError("rho not full rank in this check")
    return float(np.real(np.trace(r @ (lr - logm_h(s)))))


def ptrace(M, a, b, keep):
    T = M.reshape(a, b, a, b)
    return np.einsum('ijkj->ik', T) if keep == 'A' else np.einsum('ijil->jl', T)


def H(p, q):
    return float(sum(pi * (math.log(pi) - math.log(qi)) for pi, qi in zip(p, q) if pi > 0))


def sup_gap(rho, sig, a, b):
    lr = np.sort(np.linalg.eigvalsh(rho))[::-1]; ls = np.sort(np.linalg.eigvalsh(sig))
    rA, rB = ptrace(rho, a, b, 'A'), ptrace(rho, a, b, 'B'); sA, sB = ptrace(sig, a, b, 'A'), ptrace(sig, a, b, 'B')
    dn = lambda X: np.sort(np.linalg.eigvalsh(X))[::-1]
    return H(lr, ls) - H(dn(rA), dn(sA)) - H(dn(rB), dn(sB))


def haar(n):
    Z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / math.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q @ np.diag(np.diag(R) / np.abs(np.diag(R)))


def delta_closed(a, b, t):
    return ((1 / a) * math.log(1 - b * (a - 1) * t) + (1 / b) * math.log(1 - a * (b - 1) * t)
            - (1 / (a * b)) * math.log(1 - (a * b - 1) * t) + ((a - 1) / a) * math.log(b)
            + ((b - 1) / b) * math.log(a) + ((a - 1) * (b - 1) / (a * b)) * math.log(t))


def main():
    fails = []
    # 1. all-dimensional family
    worst = 0.0
    for a, b in [(2, 2), (2, 3), (3, 2), (3, 3), (2, 5), (4, 4)]:
        d = a * b
        for t in [1 / (4 * d), 1 / (16 * d), 1e-4, 1e-8]:
            rho = np.eye(d) / d
            e0 = np.zeros(d); e0[0] = 1
            sig = t * np.eye(d) + (1 - d * t) * np.outer(e0, e0)
            direct = D(rho, sig) - D(ptrace(rho, a, b, 'A'), ptrace(sig, a, b, 'A')) - D(ptrace(rho, a, b, 'B'), ptrace(sig, a, b, 'B'))
            closed = delta_closed(a, b, t)
            worst = max(worst, abs(direct - closed))
            if t == 1e-8 and closed >= 0:
                fails.append(('family-not-negative', a, b))
        print(f"a={a} b={b}: Delta(1e-4)={delta_closed(a, b, 1e-4):.4f}  Delta(1e-8)={delta_closed(a, b, 1e-8):.4f}")
    print(f"closed form vs direct matrix evaluation: max |error| = {worst:.2e}")
    if worst > 1e-9: fails.append(('closed-form', worst))
    t = 1 / 64
    lhs = 16 * t * (1 - 2 * t) ** 4 / (1 - 3 * t)
    print(f"a=b=2, t=1/64: exp(4 Delta) = {lhs:.6f} (< 1), Delta = {delta_closed(2, 2, t):.6f}")
    if not lhs < 1: fails.append('t=1/64')
    # 2. nondegenerate two-qubit family
    Hm = np.array([[0, 0, 0, 3], [0, 1, 0, 0], [0, 0, -1, 0], [3, 0, 0, 0]], float)
    best = -math.inf; mins = []
    for eps in [1e-6, 1e-4, 1 / 2048, 1 / 1024]:
        for eta in [1e-7, 1e-5, 1 / 8192, 1 / 4096]:
            rho = np.eye(4) / 4 + eps * Hm
            sig = np.diag([61 / 64, 1 / 64 - eta, 1 / 64, 1 / 64 + eta])
            g = sup_gap(rho, sig, 2, 2)
            best = max(best, g)
            ev = np.linalg.eigvalsh(rho); evs = np.linalg.eigvalsh(sig)
            simple = min(np.diff(np.sort(ev))) > 0 and min(np.diff(np.sort(evs))) > 0
            comm = np.linalg.norm(rho @ sig - sig @ rho)
            if not (simple and comm > 0): fails.append(('nondegenerate', eps, eta))
            # random unitaries never exceed the spectral supremum
            for _ in range(20):
                U, V, W = haar(4), haar(2), haar(2)
                gr = D(U @ rho @ U.conj().T, sig) - D(V @ ptrace(rho, 2, 2, 'A') @ V.conj().T, ptrace(sig, 2, 2, 'A')) \
                    - D(W @ ptrace(rho, 2, 2, 'B') @ W.conj().T, ptrace(sig, 2, 2, 'B'))
                if gr > g + 1e-12: fails.append(('random-exceeds-sup', eps, eta))
    print(f"nondegenerate family: largest sup_U G over the grid = {best:.6f} (claimed < -1/256 = {-1/256:.6f})")
    if not best < -1 / 256: fails.append(('nondegenerate-gap', best))
    # 3. control: product reference gives sup >= 0
    for a, b in [(2, 2), (2, 3)]:
        sA = np.diag(rng.dirichlet(np.ones(a))); sB = np.diag(rng.dirichlet(np.ones(b)))
        g = sup_gap(np.eye(a * b) / (a * b), np.kron(sA, sB), a, b)
        print(f"control product sigma a={a} b={b}: sup G = {g:.3e} (expected >= 0)")
        if g < -1e-12: fails.append(('control', a, b))
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:4]}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
