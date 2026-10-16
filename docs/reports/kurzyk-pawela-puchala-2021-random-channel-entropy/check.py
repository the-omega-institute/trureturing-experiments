"""Numerical check of Kurzyk-Pawela-Puchala, Entropy 23 (2021) 1028, Conjecture 1.

For a channel Phi on M_d with unnormalized Choi matrix C (reference first,
Tr_output C = I), put X = d*C, f(x) = x log x, A = (1/d) Tr_output f(X) and
m = Tr(A)/d.  Deterministic sandwich (checked here on samples):

    m = D at the maximally entangled input <= D(Phi || R_d) <= lambda_max(A),
    m = 2 log d - S(J),   J = C/d,
    D(omega_rho || rho x I/d) <= Tr(rho A) for every input rho.

Conjecture 1 in the additive reading: for the uniform (Hilbert-Schmidt,
environment d^2) ensemble, D(Phi || R_d) -> 1/2, i.e. H - log d -> -1/2 and
H^K - log d - H -> 0.  Prediction registered before the run: m and
lambda_max(A) both approach 1/2 and lambda_max(A) - m decreases to 0.

Usage: python3 check.py [seed]   (numpy only)
"""
import sys
import numpy as np


def sample_choi(d, rng):
    n = d * d
    g = (rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))) / np.sqrt(2)
    w = g @ g.conj().T
    t = np.einsum('iaja->ij', w.reshape(d, d, d, d))
    ev, u = np.linalg.eigh(t)
    t_mhalf = (u * ev ** -0.5) @ u.conj().T
    k = np.kron(t_mhalf, np.eye(d))
    c = k @ w @ k
    return (c + c.conj().T) / 2


def ptrace_out(m, d):
    return np.einsum('iaja->ij', m.reshape(d, d, d, d))


def fun(m, g):
    ev, u = np.linalg.eigh(m)
    return (u * g(np.clip(ev, 0, None))) @ u.conj().T


def xlogx(x):
    return np.where(x > 0, x * np.log(np.where(x > 0, x, 1)), 0.0)


def entropy(rho):
    ev = np.clip(np.linalg.eigvalsh(rho), 0, None)
    ev = ev[ev > 1e-300]
    return float(-(ev * np.log(ev)).sum())


def readings(c, d):
    x = d * c
    a = ptrace_out(fun(x, xlogx), d) / d
    a = (a + a.conj().T) / 2
    ev = np.linalg.eigvalsh(a)
    m = float(np.trace(a).real) / d
    s_j = entropy(c / d)
    return a, m, float(ev[-1]), float(ev[0]), s_j


def divergence(c, d, rho):
    ev, u = np.linalg.eigh(rho)
    sq = (u * np.sqrt(np.clip(ev, 0, None))) @ u.conj().T
    k = np.kron(sq, np.eye(d))
    omega = k @ c @ k
    return np.log(d) + entropy(rho) - entropy((omega + omega.conj().T) / 2)


def rho_of(b):
    r = b @ b.conj().T
    return r / np.trace(r).real


def optimize(c, d, rng, restarts=6, steps=600):
    best = divergence(c, d, np.eye(d) / d)
    for _ in range(restarts):
        b = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
        cur = divergence(c, d, rho_of(b))
        scale = 0.5
        for _ in range(steps):
            nb = b + scale * (rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d)))
            val = divergence(c, d, rho_of(nb))
            if val > cur:
                b, cur = nb, val
            else:
                scale *= 0.995
        best = max(best, cur)
    return best


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261016
    rng = np.random.default_rng(seed)
    ok = True
    print(f"seed {seed}")
    print("== deterministic sandwich on samples (optimized D over inputs) ==")
    for d in (2, 3, 4):
        for s in range(3):
            c = sample_choi(d, rng)
            assert np.allclose(ptrace_out(c, d), np.eye(d), atol=1e-9)
            a, m, lmax, lmin, s_j = readings(c, d)
            ident = abs(m - (2 * np.log(d) - s_j))
            d_me = divergence(c, d, np.eye(d) / d)
            d_opt = optimize(c, d, rng)
            # Tr(rho A) bound on random inputs
            worst = -1e9
            for _ in range(200):
                b = rng.standard_normal((d, d)) + 1j * rng.standard_normal((d, d))
                r = rho_of(b)
                worst = max(worst, divergence(c, d, r) - float(np.trace(r @ a).real))
            good = ident < 1e-9 and abs(d_me - m) < 1e-9 and m - 1e-9 <= d_opt <= lmax + 1e-9 and worst <= 1e-9
            ok &= good
            print(f"d={d} s={s} m={m:.6f} D_opt={d_opt:.6f} lmax(A)={lmax:.6f} "
                  f"|m-(2logd-S(J))|={ident:.1e} max[D-Tr(rhoA)]={worst:.2e} {'ok' if good else 'FAIL'}")
    print("== asymptotics: m, lambda_min(A), lambda_max(A) ==")
    for d in (2, 4, 6, 8, 12, 16, 24, 32):
        reps = 20 if d <= 12 else (6 if d <= 24 else 3)
        rows = np.array([readings(sample_choi(d, rng), d)[1:4] for _ in range(reps)])
        mm, lmax, lmin = rows.mean(axis=0)
        print(f"d={d:2d} reps={reps:2d} mean m={mm:.5f} mean lmin(A)={lmin:.5f} "
              f"mean lmax(A)={lmax:.5f} mean gap bound lmax-m={(rows[:,1]-rows[:,0]).mean():.5f}")
    print("SANDWICH_OK" if ok else "SANDWICH_FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
