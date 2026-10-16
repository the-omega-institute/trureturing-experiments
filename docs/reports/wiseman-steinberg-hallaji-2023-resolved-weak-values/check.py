"""Numerical check of the Wiseman-Steinberg-Hallaji conjecture (arXiv:2211.15761v2, Sec. V).

Passive linear optics on M modes: input coherent state |alpha, 0, ...>, network A,
weak measurement of the photon number n_1 of mode 1, network B, photon-number
resolving detection of m photons in mode 1 (other modes ignored).  With
E_m = |m><m| (x) I,
    W_m = <B A psi| E_m B n_1 A psi> / <B A psi| E_m |B A psi>      (source Eq. (2)),
and w = the same weak value for a single-photon input post-selected on a click (m = 1).
Conjecture: W_m - W_0 = m * w for every m >= 1 with nonzero probability.

The Fock representation Gamma(U) is built independently of any closed formula, by
exponentiating the second-quantized generator a^dagger H a in a truncated Fock space.
Prediction registered before the run: |W_m - W_0 - m w| is at truncation level
(< 1e-6) for random networks with M = 2, 3, and m = 1..4.

Usage: python3 check.py [seed [M,N ...]]   (numpy; default truncations 2,18 and 3,11)
"""
import sys
import itertools
import numpy as np


def fock_ops(M, N):
    states = list(itertools.product(range(N + 1), repeat=M))
    idx = {s: i for i, s in enumerate(states)}
    D = len(states)
    a = []
    for j in range(M):
        A = np.zeros((D, D), complex)
        for s in states:
            if s[j] > 0:
                t = list(s); t[j] -= 1
                A[idx[tuple(t)], idx[s]] = np.sqrt(s[j])
        a.append(A)
    return states, idx, a


def gamma(U, a):
    """Gamma(U) = exp(i a^dag H a) with U = exp(i H), restricted to the truncated space."""
    ev, V = np.linalg.eig(U)
    H = V @ np.diag(np.angle(ev)) @ np.linalg.inv(V)
    H = (H + H.conj().T) / 2
    G = sum(H[j, k] * a[j].conj().T @ a[k] for j in range(len(a)) for k in range(len(a)))
    w, Q = np.linalg.eigh(G)
    return Q @ np.diag(np.exp(1j * w)) @ Q.conj().T


def random_unitary(M, rng):
    z = (rng.standard_normal((M, M)) + 1j * rng.standard_normal((M, M))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    return q * (np.diag(r) / np.abs(np.diag(r)))


def coherent_vec(beta, states):
    from math import factorial, sqrt
    v = np.array([np.prod([b ** n / sqrt(factorial(n)) for b, n in zip(beta, s)]) for s in states], complex)
    return v * np.exp(-np.sum(np.abs(beta) ** 2) / 2)


def projector(m, states):
    return np.diag([1.0 if s[0] == m else 0.0 for s in states]).astype(complex)


def weak_value(GA, GB, E, nop, psi):
    psi_p = GA @ psi
    fin = GB @ psi_p
    num = np.vdot(fin, E @ (GB @ (nop @ psi_p)))
    den = np.vdot(fin, E @ fin)
    return num / den, den.real


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261016
    rng = np.random.default_rng(seed)
    ok = True
    worst = 0.0
    TRUNC = ((2, 18), (3, 11)) if len(sys.argv) <= 2 else tuple(tuple(int(x) for x in s.split(",")) for s in sys.argv[2:])
    for M, N in TRUNC:
        states, idx, a = fock_ops(M, N)
        nop = a[0].conj().T @ a[0]
        single = np.zeros(len(states), complex)
        single[idx[(1,) + (0,) * (M - 1)]] = 1
        for trial in range(4):
            A = random_unitary(M, rng); B = random_unitary(M, rng)
            GA, GB = gamma(A, a), gamma(B, a)
            alpha = complex(rng.uniform(0.4, 1.0), rng.uniform(-0.6, 0.6))
            psi = coherent_vec(np.array([alpha] + [0] * (M - 1)), states)
            w, pw = weak_value(GA, GB, projector(1, states), nop, single)
            W0, p0 = weak_value(GA, GB, projector(0, states), nop, psi)
            line = f"M={M} trial={trial} alpha={alpha:.3f} w={w:.6f} (click prob {pw:.4f}) W0={W0:.6f}"
            print(line)
            for m in range(1, 5):
                Wm, pm = weak_value(GA, GB, projector(m, states), nop, psi)
                err = abs(Wm - W0 - m * w)
                worst = max(worst, err)
                ok &= err < 1e-6
                print(f"   m={m}: P={pm:.3e}  W_m - W_0 = {Wm - W0:.8f}  m*w = {m * w:.8f}  |diff| = {err:.2e}")
    print(f"largest |W_m - W_0 - m w| = {worst:.2e}")
    print("ALL_OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
