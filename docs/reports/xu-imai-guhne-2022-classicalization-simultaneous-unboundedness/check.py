"""Xu-Imai-Guhne, "Fate of multiparticle entanglement when one particle becomes classical", arXiv:2206.12834v2,
Section V.B question: can Delta_E and tildeDelta_E - Delta_E both be arbitrarily large with C a qubit?

Search-seat answer: E = N_ABC (Eq. (4)), rho_d = Phi_d (x) W on A = A'A0, B = B'B0, C (dim A' = dim B' = d).
Claims: E(rho_d) = d - 1 + (2d+1) sqrt(2)/3; sup over entanglement-breaking channels on C of E(Phi(rho_d)) = 5d/3 - 1,
attained by the computational-basis measurement with orthogonal register states; hence
Delta = (2 sqrt 2 - 2) d / 3 + sqrt 2 / 3 and tildeDelta - Delta = 5d/3 - 1 - E(rho_AB (x) |0><0|), both linear in d.

Checks (floating point, numpy, tolerance 1e-9):
 1. E(rho_d) for d = 2, 3, 4 against the formula.
 2. The attaining channel: E(Phi(rho_d)) = 5d/3 - 1.
 3. 400 random channels per d (random POVMs with 2-4 rank-one or full-rank effects, random register states of
    dimension 1-3): E(Phi(rho_d)) <= 5d/3 - 1 + 1e-9.
 4. tildeDelta - Delta > 0 and both grow with d (values printed).
"""
import numpy as np
rng = np.random.default_rng(220612834)

def ptrans(rho, dims, sub):
    # partial transpose on subsystem index `sub` of a multipartite matrix with local dims `dims`
    n = len(dims)
    r = rho.reshape(dims + dims)
    perm = list(range(2 * n))
    perm[sub], perm[n + sub] = perm[n + sub], perm[sub]
    return r.transpose(perm).reshape(rho.shape)

def tnorm(M):
    return np.linalg.svd(M, compute_uv=False).sum()

def neg(rho, dims, sub):
    return (tnorm(ptrans(rho, dims, sub)) - 1) / 2

def NABC(rho, dims3):
    # dims3 = (dA, dB, dC); N_{AB|C} + N_{BC|A} + N_{AC|B} = transposes on C, A, B
    return sum(neg(rho, list(dims3), s) for s in (2, 0, 1))

def state(d):
    phi = np.zeros(d * d); 
    for r in range(d): phi[r * d + r] = 1
    phi /= np.sqrt(d)
    W = np.zeros(8); W[1] = W[2] = W[4] = 1; W /= np.sqrt(3)   # |001>,|010>,|100> on A0 B0 C
    # order factors as A' B' A0 B0 C, then regroup A = A'A0, B = B'B0
    psi = np.kron(phi, W)  # A' B' A0 B0 C
    t = psi.reshape(d, d, 2, 2, 2).transpose(0, 2, 1, 3, 4).reshape(-1)  # A' A0 B' B0 C
    return np.outer(t, t.conj())

def classicalize(rho, dA, dB, effects, taus):
    dE = taus[0].shape[0]
    r = rho.reshape(dA, dB, 2, dA, dB, 2)
    out = np.zeros((dA * dB * dE, dA * dB * dE), dtype=complex)
    for m, tau in zip(effects, taus):
        br = np.einsum('abcxyz,zc->abxy', r, m).reshape(dA * dB, dA * dB)
        out += np.kron(br, tau)
    return out, dE

def rand_povm(k):
    G = [rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)) for _ in range(k)]
    P = [g @ g.conj().T for g in G]
    S = sum(P)
    w, V = np.linalg.eigh(S); Sm = V @ np.diag(w ** -0.5) @ V.conj().T
    return [Sm @ p @ Sm for p in P]

def rand_state(e):
    g = rng.normal(size=(e, e)) + 1j * rng.normal(size=(e, e)); s = g @ g.conj().T; return s / np.trace(s)

ok = True
for d in (2, 3, 4):
    rho = state(d); dA = dB = 2 * d
    E0 = NABC(rho, (dA, dB, 2))
    f = d - 1 + (2 * d + 1) * np.sqrt(2) / 3
    sup = 5 * d / 3 - 1
    eff = [np.diag([1, 0]).astype(complex), np.diag([0, 1]).astype(complex)]
    taus = [np.diag([1, 0]).astype(complex), np.diag([0, 1]).astype(complex)]
    out, dE = classicalize(rho, dA, dB, eff, taus)
    att = NABC(out, (dA, dB, dE))
    mx = -1
    for t in range(400):
        k = rng.integers(2, 5); e = rng.integers(1, 4)
        P = rand_povm(k); T = [rand_state(e) for _ in range(k)]
        o, de = classicalize(rho, dA, dB, P, T)
        mx = max(mx, NABC(o, (dA, dB, de)))
    rhoAB = np.einsum('abcxyc->abxy', rho.reshape(dA, dB, 2, dA, dB, 2)).reshape(dA * dB, dA * dB)
    tilde_out = np.kron(rhoAB, np.diag([1, 0]).astype(complex))
    Et = NABC(tilde_out, (dA, dB, 2))
    Delta = E0 - sup; tD = E0 - Et
    print(f"d={d}: E(rho)={E0:.9f} formula={f:.9f}  attained={att:.9f} 5d/3-1={sup:.9f}  random max={mx:.9f}  "
          f"Delta={Delta:.6f} tildeDelta-Delta={tD - Delta:.6f}")
    ok &= abs(E0 - f) < 1e-9 and abs(att - sup) < 1e-9 and mx <= sup + 1e-9 and Delta > 0 and tD - Delta > 0
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
