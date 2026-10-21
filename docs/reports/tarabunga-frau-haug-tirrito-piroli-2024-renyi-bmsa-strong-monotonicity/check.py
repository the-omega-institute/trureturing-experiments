"""Rényi BMSA strong monotonicity checks (Tarabunga-Frau-Haug-Tirrito-Piroli, arXiv:2411.05766, Eqs. (5), (31)).
A_alpha(psi) = min over stabilizer bases of the Rényi entropy of the measurement distribution of psi in that basis
(the dephasing G_{G(s)}(rho) of Eq. (31)). Strong monotonicity (5): A(psi) >= sum_lambda p_lambda A(psi_lambda) for
computational-basis measurements of a proper qubit subset, outputs normalized and keeping the measured qubits."""
import itertools, sys
import numpy as np

I2 = np.eye(2); X = np.array([[0, 1], [1, 0]]); Z = np.diag([1, -1])

def pauli(x, z):
    m = np.array([[1.0 + 0j]])
    for xi, zi in zip(x, z):
        p = np.eye(2, dtype=complex)
        if xi: p = p @ X
        if zi: p = p @ Z
        m = np.kron(m, p)
    ph = (1j) ** (sum(a & b for a, b in zip(x, z)) % 4)
    return ph * m  # Hermitian: i^{x.z} X^x Z^z

def lagrangians(N):
    out = []
    vecs = [v for v in itertools.product([0, 1], repeat=2 * N) if any(v)]
    def symp(u, v):
        return sum(u[i] * v[N + i] + u[N + i] * v[i] for i in range(N)) % 2
    def rref_key(rows):
        M = [list(r) for r in rows]
        r = 0
        for c in range(2 * N):
            piv = next((i for i in range(r, len(M)) if M[i][c]), None)
            if piv is None: continue
            M[r], M[piv] = M[piv], M[r]
            for i in range(len(M)):
                if i != r and M[i][c]:
                    M[i] = [(a + b) % 2 for a, b in zip(M[i], M[r])]
            r += 1
        return tuple(tuple(x) for x in M[:r])
    seen = set()
    def extend(rows):
        if len(rows) == N:
            k = rref_key(rows)
            if k not in seen:
                seen.add(k); out.append(k)
            return
        span = set()
        for coeffs in itertools.product([0, 1], repeat=len(rows)):
            s = tuple(sum(c * r[i] for c, r in zip(coeffs, rows)) % 2 for i in range(2 * N))
            span.add(s)
        start = vecs.index(rows[-1]) + 1 if rows else 0
        for v in vecs[start:]:
            if v in span: continue
            if all(symp(v, r) == 0 for r in rows):
                extend(rows + [v])
    extend([])
    return out

def distributions(N, L, psi):
    gens = [pauli(r[:N], r[N:]) for r in L]
    # signs so generators commute (they do); product over subsets
    d = 2 ** N
    exp = {}
    for T in itertools.product([0, 1], repeat=N):
        M = np.eye(d, dtype=complex)
        for t, g in zip(T, gens):
            if t: M = M @ g
        exp[T] = np.vdot(psi, M @ psi)
    p = []
    for s in itertools.product([1, -1], repeat=N):
        val = 0
        for T, e in exp.items():
            sign = np.prod([si if t else 1 for si, t in zip(s, T)])
            val += sign * e
        p.append((val / d).real)
    return np.clip(np.array(p), 0, None)

def renyi(p, a):
    p = p[p > 1e-14]
    if abs(a - 1) < 1e-12: return float(-(p * np.log(p)).sum())
    if a == 0: return float(np.log(len(p)))
    return float(np.log((p ** a).sum()) / (1 - a))

def bmsa(N, Ls, psi, alphas):
    best = {a: np.inf for a in alphas}
    for L in Ls:
        p = distributions(N, L, psi)
        for a in alphas:
            best[a] = min(best[a], renyi(p, a))
    return best

def measure(N, psi, subset):
    outs = []
    for lam in itertools.product([0, 1], repeat=len(subset)):
        phi = psi.copy().reshape([2] * N)
        idx = [slice(None)] * N
        mask = np.zeros([2] * N, dtype=bool)
        sl = [slice(None)] * N
        for q, b in zip(subset, lam): sl[q] = b
        mask[tuple(sl)] = True
        phi = np.where(mask, phi, 0).reshape(-1)
        pr = np.vdot(phi, phi).real
        if pr > 1e-12: outs.append((pr, phi / np.sqrt(pr)))
    return outs

def main():
    rng = np.random.default_rng(20261021)
    alphas = [0.0, 0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0]
    worst = {a: (-np.inf, None) for a in alphas}
    for N, trials in [(2, 300), (3, 120)]:
        Ls = lagrangians(N)
        print(f"N={N}: stabilizer bases {len(Ls)}")
        subsets = [s for k in range(1, N) for s in itertools.combinations(range(N), k)]
        for t in range(trials):
            kind = t % 3
            if kind == 0:
                psi = rng.normal(size=2 ** N) + 1j * rng.normal(size=2 ** N)
            elif kind == 1:  # sparse, near-stabilizer superpositions
                psi = np.zeros(2 ** N, complex); idx = rng.choice(2 ** N, size=2, replace=False)
                psi[idx] = rng.normal(size=2) + 1j * rng.normal(size=2)
            else:
                th = rng.uniform(0, np.pi, size=N)
                v = np.array([1.0 + 0j])
                for x in th: v = np.kron(v, np.array([np.cos(x / 2), np.sin(x / 2)]))
                psi = v + 0.3 * (rng.normal(size=2 ** N) + 1j * rng.normal(size=2 ** N))
            psi = psi / np.linalg.norm(psi)
            A = bmsa(N, Ls, psi, alphas)
            for sub in subsets:
                outs = measure(N, psi, sub)
                avg = {a: 0.0 for a in alphas}
                for pr, phi in outs:
                    B = bmsa(N, Ls, phi, alphas)
                    for a in alphas: avg[a] += pr * B[a]
                for a in alphas:
                    gap = avg[a] - A[a]
                    if gap > worst[a][0]: worst[a] = (gap, (N, t, sub))
    bad = 0
    for a in alphas:
        g, where = worst[a]
        verdict = "violation" if g > 1e-9 else "holds"
        if a <= 1 and g > 1e-9: bad += 1
        print(f"alpha={a}: max(avg_out - A_in) = {g:.3e} ({verdict}) at {where}")
    # reproduction of the source's numbered bounds: (33) 0 <= A <= N ln 2, faithfulness, (51) D_min <= A_alpha,
    # (43) M_n <= n/(n-1) A_alpha for n > 1, alpha <= 2
    N = 3; Ls = lagrangians(N); d = 2 ** N
    allp = [pauli(x, z) for x in itertools.product([0, 1], repeat=N) for z in itertools.product([0, 1], repeat=N)]
    stabvecs = []
    for L in Ls:
        gens = [pauli(r[:N], r[N:]) for r in L]
        for s in itertools.product([1, -1], repeat=N):
            P = np.eye(d, dtype=complex)
            for si, g in zip(s, gens): P = P @ (np.eye(d) + si * g) / 2
            w, v = np.linalg.eigh(P)
            stabvecs.append(v[:, -1])
    rep_bad = 0
    for t in range(60):
        psi = rng.normal(size=d) + 1j * rng.normal(size=d); psi /= np.linalg.norm(psi)
        A = bmsa(N, Ls, psi, [0.5, 1.0, 1.5, 2.0])
        dmin = -np.log(max(abs(np.vdot(s, psi)) ** 2 for s in stabvecs))
        r = np.array([np.vdot(psi, P @ psi).real for P in allp])
        for a in [0.5, 1.0, 1.5, 2.0]:
            rep_bad += not (-1e-9 <= A[a] <= N * np.log(2) + 1e-9)
            rep_bad += not (dmin <= A[a] + 1e-9)
        for nn in [1.5, 2.0, 3.0]:
            M = np.log((np.abs(r) ** (2 * nn)).sum() / d) / (1 - nn)
            for a in [0.5, 1.0, 1.5, 2.0]:
                rep_bad += not (M <= nn / (nn - 1) * A[a] + 1e-9)
    zero = max(max(bmsa(N, Ls, s, [0.5, 1.0, 2.0]).values()) for s in stabvecs[:40])
    print(f"source bounds (33), (51), (43) on 60 random 3-qubit states: violations {rep_bad}; max A over 40 stabilizer states {zero:.2e}")
    bad += rep_bad + (zero > 1e-9)
    # negative family arithmetic (alpha = 3/2, n = 29)
    c2 = (2 + np.sqrt(2)) / 4; s2 = (2 - np.sqrt(2)) / 4; a = 1.5; n = 29
    up = np.log(2 ** (-a) * (1 + (c2 ** a + s2 ** a) ** n)) / (1 - a)
    low = n / 4 * np.log(4 / 3)
    print(f"family alpha=3/2 n=29: A_in <= {up:.6f}, avg_out >= {low:.6f}, 3 log 2 = {3*np.log(2):.6f}, 2^46 > 3^29: {2**46 > 3**29}")
    bad += not (up < 3 * np.log(2) < low)
    print(f"EXIT={1 if bad else 0}")
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
