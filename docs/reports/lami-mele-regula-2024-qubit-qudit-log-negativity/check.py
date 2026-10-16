"""Numerical check for Lami-Mele-Regula, arXiv:2405.09613v2, Remark S29.

On a 2 x D system the exact PPT entanglement cost equals E_kappa = log2 kappa
(their Corollary 2), with kappa(rho) = min{Tr S : S >= 0, -S^G <= rho^G <= S^G}.
Remark S29 asks whether Ishizaka's two-qubit theorem (|rho^G|^G >= 0, hence
E_kappa = E_N) extends to all qubit-qudit systems.

Family: u = a|00> + b|01>, v = c|02> + d|10>, rho = (uu* + vv*)/Z.
Prediction registered before the run: for every a, b, c, d > 0 the matrix
|rho^G|^G has a negative eigenvalue, and kappa(rho) > Tr|rho^G| strictly
(primal and dual SDP values agree and exceed the trace norm); the same holds
after mixing with white noise in 2 x D for D = 3, 4, 6.

Usage: python3 check.py [seed]   (numpy, cvxpy with the CLARABEL solver)
"""
import sys
import numpy as np
import cvxpy as cp


def pt(m, da, db):
    return m.reshape(da, db, da, db).transpose(0, 3, 2, 1).reshape(da * db, da * db)


def absm(h):
    ev, u = np.linalg.eigh(h)
    return (u * np.abs(ev)) @ u.conj().T


def family(a, b, c, d, D=3):
    def k(i, j):
        x = np.zeros(2 * D)
        x[i * D + j] = 1
        return x
    u = a * k(0, 0) + b * k(0, 1)
    v = c * k(0, 2) + d * k(1, 0)
    r = np.outer(u, u) + np.outer(v, v)
    return r / np.trace(r)


def pt_expr(x, da, db):
    # partial transpose of a cvxpy matrix expression on the second factor
    n = da * db
    blocks = [[None] * da for _ in range(da)]
    for i in range(da):
        for k in range(da):
            blocks[i][k] = x[i * db:(i + 1) * db, k * db:(k + 1) * db].T
    return cp.bmat(blocks)


def kappa(rho, D):
    n = 2 * D
    h = pt(rho, 2, D)
    s = cp.Variable((n, n), symmetric=True)
    sg = pt_expr(s, 2, D)
    cons = [s >> 0, sg - h >> 0, sg + h >> 0]
    primal = cp.Problem(cp.Minimize(cp.trace(s)), cons)
    primal.solve(solver=cp.CLARABEL)
    p = cp.Variable((n, n), symmetric=True)
    q = cp.Variable((n, n), symmetric=True)
    dual = cp.Problem(cp.Maximize(cp.trace(h @ (p - q))),
                      [p >> 0, q >> 0, np.eye(n) - pt_expr(p + q, 2, D) >> 0])
    dual.solve(solver=cp.CLARABEL)
    return primal.value, dual.value


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261016
    rng = np.random.default_rng(seed)
    ok = True
    print(f"seed {seed}")
    rho = family(np.sqrt(2), 1, 1, np.sqrt(2))
    h = pt(rho, 2, 3)
    print("seed a=d=sqrt2, b=c=1: 6*eig(rho^G) =", np.round(6 * np.linalg.eigvalsh(h), 6))
    b = pt(absm(h), 2, 3)
    tn = np.abs(np.linalg.eigvalsh(h)).sum()
    kp, kd = kappa(rho, 3)
    print(f"  min eig |rho^G|^G = {np.linalg.eigvalsh(b)[0]:.6f}; Tr|rho^G| = {tn:.6f} (4/3 = {4/3:.6f}); "
          f"kappa primal {kp:.6f} dual {kd:.6f}; E_N = {np.log2(tn):.6f} E_kappa = {np.log2(kp):.6f}")
    ok &= np.linalg.eigvalsh(b)[0] < -1e-6 and kd > tn + 1e-6 and abs(kp - kd) < 1e-5
    print("== random family members (D = 3) ==")
    worst_gap = np.inf
    for t in range(40):
        a, bb, c, d = rng.uniform(0.05, 3, 4)
        rho = family(a, bb, c, d)
        h = pt(rho, 2, 3)
        mb = np.linalg.eigvalsh(pt(absm(h), 2, 3))[0]
        tn = np.abs(np.linalg.eigvalsh(h)).sum()
        kp, kd = kappa(rho, 3)
        gap = kd - tn
        worst_gap = min(worst_gap, gap)
        ok &= mb < 0 and gap > 0 and abs(kp - kd) < 1e-5
    print(f"  40 samples: every binegativity eigenvalue < 0 and kappa > Tr|rho^G|; smallest dual gap {worst_gap:.2e}")
    print("== white noise, 2 x D ==")
    for D in (3, 4, 6):
        for eps in (0.01, 0.1):
            r0 = family(np.sqrt(2), 1, 1, np.sqrt(2), D)
            rho = (1 - eps) * r0 + eps * np.eye(2 * D) / (2 * D)
            h = pt(rho, 2, D)
            mb = np.linalg.eigvalsh(pt(absm(h), 2, D))[0]
            tn = np.abs(np.linalg.eigvalsh(h)).sum()
            kp, kd = kappa(rho, D)
            good = mb < 0 and kd > tn + 1e-7 and abs(kp - kd) < 1e-5
            ok &= good
            print(f"  D={D} eps={eps}: min eig rho {np.linalg.eigvalsh(rho)[0]:.5f}, min eig |rho^G|^G {mb:.6f}, "
                  f"Tr|rho^G| {tn:.6f}, kappa {kp:.6f}/{kd:.6f} {'ok' if good else 'FAIL'}")
    print("== control: two qubits (Ishizaka) ==")
    for t in range(20):
        g = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
        rho = g @ g.conj().T
        rho /= np.trace(rho).real
        rho = rho.real if False else rho
        h = pt(rho, 2, 2)
        mb = np.linalg.eigvalsh(pt(absm(h), 2, 2))[0]
        ok &= mb > -1e-9
    print("  20 random two-qubit states: binegativity operator PSD")
    print("ALL_OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
