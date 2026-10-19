"""Badhani-Das, arXiv:2604.01217v2, Section VI.1 conjecture (after Proposition 10).

Claim checked: for N_p = p R^pi + (1 - p) U with U a unitary on C^d (x) C^2 and B' = B = C^2,
S_inf[A|B]_{N_p} = S_inf(R_A A | R_B B)_{Phi^{N_p}} - log|A'|.

Both sides are computed by independent SDPs (cvxpy, SCS and CLARABEL):
  state side   q_state = min{Tr M : omega <= 1_X (x) M, M >= 0},       S_inf(X|Y) = -log2 q_state;
  channel side q_chan  = min{Tr M / 2 : Gamma <= 1_X (x) M, M >= 0, Tr_B M = (Tr M / 2) 1_{R_B}},
               S_inf[A|B] = -log2 q_chan   (source Eq. con-min-ent-sdp).
The identity is q_chan = d * q_state.  Negative controls: random CPTP maps and unitaries with
non-white (pure-state replacement) noise, where the identity is expected to fail in general.
Also checks the structural steps of the proof: a maximally entangled eigenbasis of rho_Y and
feasibility of the pinched optimizer.

Requires numpy, scipy, cvxpy.  Usage: python3 check.py  (prints RESULT PASS / RESULT FAIL).
"""
import numpy as np
import cvxpy as cp
from scipy.stats import unitary_group

rng = np.random.default_rng(20261019)
TOL = 2e-5


def choi_state(kraus_or_fn, d):
    """Normalised Choi state on (R_A A)(R_B B), index order (ra, a, rb, b)."""
    n = 2 * d
    phi = np.zeros((n * n, n * n), dtype=complex)
    # input basis |a', b'>, index a'*2+b'; reference copies r_a, r_b
    for i in range(n):
        for j in range(n):
            Eij = np.zeros((n, n), dtype=complex); Eij[i, j] = 1
            out = kraus_or_fn(Eij)  # on (a, b)
            ra, rb = divmod(i, 2); sa, sb = divmod(j, 2)
            for a in range(d):
                for b in range(2):
                    for a2 in range(d):
                        for b2 in range(2):
                            row = ((ra * d + a) * 2 + rb) * 2 + b
                            col = ((sa * d + a2) * 2 + sb) * 2 + b2
                            phi[row, col] += out[a * 2 + b, a2 * 2 + b2]
    return phi / n


def state_sdp(omega, dX):
    M = cp.Variable((4, 4), hermitian=True)
    cons = [M >> 0, cp.kron(np.eye(dX), M) - omega >> 0]
    prob = cp.Problem(cp.Minimize(cp.real(cp.trace(M))), cons)
    return solve(prob)


def channel_sdp(Gamma, dX):
    M = cp.Variable((4, 4), hermitian=True)
    t = cp.real(cp.trace(M)) / 2
    trB = cp.partial_trace(M, [2, 2], axis=1)
    cons = [M >> 0, cp.kron(np.eye(dX), M) - Gamma >> 0, trB == t * np.eye(2)]
    prob = cp.Problem(cp.Minimize(t), cons)
    return solve(prob)


def solve(prob):
    vals = []
    for s in ("CLARABEL", "SCS"):
        try:
            kw = {"eps": 1e-9, "max_iters": 200000} if s == "SCS" else {}
            prob.solve(solver=s, **kw)
            vals.append(prob.value)
        except Exception:
            pass
    return vals[0], vals


def noisy_unitary(U, p, d, sigma=None):
    n = 2 * d
    sig = np.eye(n) / n if sigma is None else sigma
    return lambda Z: (1 - p) * U @ Z @ U.conj().T + p * np.trace(Z) * sig


def random_channel(d, k=3):
    n = 2 * d
    G = rng.normal(size=(k * n, n)) + 1j * rng.normal(size=(k * n, n))
    Q, _ = np.linalg.qr(G)
    Ks = [Q[i * n:(i + 1) * n, :] for i in range(k)]
    return lambda Z: sum(K @ Z @ K.conj().T for K in Ks)


def structural(U, d):
    """Maximally entangled eigenbasis of rho_Y and pinching feasibility (steps 1-4)."""
    omega = choi_state(noisy_unitary(U, 0.0, d), d)
    dX = d * d
    rhoY = np.einsum('xyxz->yz', omega.reshape(dX, 4, dX, 4))
    # marginals
    r = rhoY.reshape(2, 2, 2, 2)
    ok = np.allclose(np.einsum('abcb->ac', r), np.eye(2) / 2, atol=1e-10)
    ok &= np.allclose(np.einsum('abad->bd', r), np.eye(2) / 2, atol=1e-10)
    # magic basis: rho_Y is real there
    W = np.array([[1, 0, 0, 1j], [0, 1j, 1, 0], [0, 1j, -1, 0], [1, 0, 0, -1j]]) / np.sqrt(2)
    R = W.conj().T @ rhoY @ W
    ok &= np.allclose(R.imag, 0, atol=1e-10)
    w, V = np.linalg.eigh(R.real)
    F = W @ V
    for i in range(4):
        f = F[:, i].reshape(2, 2)
        ok &= np.allclose(f @ f.conj().T, np.eye(2) / 2, atol=1e-10)
    return bool(ok)


def main():
    fails = 0
    print("identity checks: q_chan vs d*q_state")
    for d in (1, 2, 3):
        for trial in range(4):
            U = unitary_group.rvs(2 * d, random_state=rng)
            for p in (0.0, 0.25, 0.6, 1.0):
                omega = choi_state(noisy_unitary(U, p, d), d)
                qs, _ = state_sdp(omega, d * d)
                qc, _ = channel_sdp(2 * d * omega, d * d)
                diff = abs(qc - d * qs)
                ok = diff < TOL * max(1, qc)
                fails += not ok
                print(f"d={d} trial={trial} p={p:.2f} q_state={qs:.8f} q_chan={qc:.8f} d*q_state={d*qs:.8f} diff={diff:.2e} {'ok' if ok else 'FAIL'}")
    # named unitaries from the source figure, d = 2
    CNOT = np.eye(4)[[0, 1, 3, 2]]; SWAP = np.eye(4)[[0, 2, 1, 3]]
    for name, U in (("CNOT", CNOT), ("SWAP", SWAP), ("id", np.eye(4))):
        for p in (0.0, 0.5, 0.9):
            omega = choi_state(noisy_unitary(U.astype(complex), p, 2), 2)
            qs, _ = state_sdp(omega, 4); qc, _ = channel_sdp(4 * omega, 4)
            diff = abs(qc - 2 * qs); ok = diff < TOL * max(1, qc); fails += not ok
            print(f"{name} p={p} S_chan={-np.log2(qc):.6f} S_state-1={-np.log2(qs)-1:.6f} {'ok' if ok else 'FAIL'}")
    print("structural steps (maximally entangled eigenbasis via magic basis):")
    for d in (1, 2, 3, 4):
        for trial in range(5):
            U = unitary_group.rvs(2 * d, random_state=rng)
            s = structural(U, d); fails += not s
            print(f"d={d} trial={trial} structural={'ok' if s else 'FAIL'}")
    print("negative controls (identity not expected in general):")
    gaps = []
    for d in (1, 2):
        for trial in range(4):
            N = random_channel(d)
            omega = choi_state(N, d)
            qs, _ = state_sdp(omega, d * d); qc, _ = channel_sdp(2 * d * omega, d * d)
            gaps.append(qc - d * qs)
            print(f"random CPTP d={d} trial={trial} q_chan-d*q_state={qc-d*qs:.6f}")
        U = unitary_group.rvs(2 * d, random_state=rng)
        sig = np.zeros((2 * d, 2 * d)); sig[0, 0] = 1
        omega = choi_state(noisy_unitary(U, 0.5, d, sig), d)
        qs, _ = state_sdp(omega, d * d); qc, _ = channel_sdp(2 * d * omega, d * d)
        gaps.append(qc - d * qs)
        print(f"pure-replacement noise d={d} q_chan-d*q_state={qc-d*qs:.6f}")
    discriminates = max(gaps) > 1e-3
    print(f"negative controls discriminate (some gap > 1e-3): {discriminates}")
    fails += not discriminates
    print("RESULT", "PASS" if fails == 0 else f"FAIL ({fails})")


if __name__ == "__main__":
    main()
