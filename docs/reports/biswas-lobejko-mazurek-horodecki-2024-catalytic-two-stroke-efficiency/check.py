"""Biswas–Łobejko–Mazurek–Horodecki arXiv:2402.10384v2 (PRE 110, 044120 (2024)), Section VII conjecture:
"catalytic enhancement is always possible" for two-stroke heat engines with d-level hot and cold systems.

For random d-level Hamiltonians and 0 < beta_h < beta_c:
1. eta_* = max over population permutations P of the 2-body state with W > 0 of W/Q_h (Birkhoff: this is the
   noncatalytic optimum over unitaries); check eta_* < eta_C = 1 - beta_h/beta_c.
2. Catalyst cycle (seat construction): levels z0=(h0,c0), z1=(h1,c0), z2=(h1,c1); k edges to z0, l to z2 with
   beta_h/beta_c < l y/(k x) < 1 - eta_*; catalyst s_i from S_i; U = product of transpositions (i,z1)<->(i+1,v_i).
   Check: U is a permutation of the joint basis, catalyst marginal restored exactly, W > 0, eta > eta_*,
   and eta = 1 - l y/(k x).
3. Random unitaries on the 2-body system never beat eta_* (sanity of the permutation reduction).
"""
import itertools, math, sys
import numpy as np

rng = np.random.default_rng(20261017)


def engine_quantities(p, q, Eh, Ec, d):
    P = p.reshape(d, d); Q = q.reshape(d, d)
    Qh = float(Eh @ (P.sum(1) - Q.sum(1))); Qc = float(Ec @ (P.sum(0) - Q.sum(0)))
    return Qh, Qc, Qh + Qc


def eta_star(p, Eh, Ec, d):
    best = -math.inf; n = d * d
    perms = np.array(list(itertools.permutations(range(n))))
    Q = p[perms]                                   # q = P p with q_j = p_{perm[j]}
    QH = (Q.reshape(-1, d, d).sum(2) @ Eh); QC = (Q.reshape(-1, d, d).sum(1) @ Ec)
    PH = p.reshape(d, d).sum(1) @ Eh; PC = p.reshape(d, d).sum(0) @ Ec
    Qh = PH - QH; Qc = PC - QC; W = Qh + Qc
    mask = W > 1e-14
    if not mask.any(): return None
    return float((W[mask] / Qh[mask]).max())


def gibbs(E, b):
    w = np.exp(-b * E); return w / w.sum()


def trial(d):
    Eh = np.sort(rng.uniform(0, 3, d)); Eh -= Eh[0]
    Ec = np.sort(rng.uniform(0, 3, d)); Ec -= Ec[0]
    bh = rng.uniform(0.1, 1.0); bc = bh * rng.uniform(1.5, 6.0)
    p = np.kron(gibbs(Eh, bh), gibbs(Ec, bc))
    es = eta_star(p, Eh, Ec, d)
    if es is None: return None
    etaC = 1 - bh / bc
    # levels: h0=0,h1=1 ; c0=0,c1=1
    x = Eh[1] - Eh[0]; y = Ec[1] - Ec[0]
    A = (bh / bc) * x / y; B = (1 - es) * x / y
    k = 1
    while k * (B - A) <= 1: k += 1
    l = math.floor(k * A) + 1
    m = k + l
    idx = lambda h, c: h * d + c
    z0, z1, z2 = idx(0, 0), idx(1, 0), idx(1, 1)
    a, b, c = p[z1], p[z0], p[z2]
    v = [z0] * k + [z2] * l; r = [b / a] * k + [c / a] * l
    R = math.prod(r)
    S = [sum(math.prod(r[(i + t) % m] for t in range(j)) for j in range(m)) for i in range(m)]
    C = sum(S); s = np.array(S) / C
    joint = np.kron(s, p)                           # index = i*d*d + z
    N = m * d * d; perm = np.arange(N)
    used = set()
    for i in range(m):
        src = i * d * d + z1; tgt = ((i + 1) % m) * d * d + v[i]
        if src in used or tgt in used: return ("not-disjoint",)
        used |= {src, tgt}; perm[src], perm[tgt] = tgt, src
    out = joint[perm]                               # permutation unitary output populations
    cat_before = joint.reshape(m, -1).sum(1); cat_after = out.reshape(m, -1).sum(1)
    body_after = out.reshape(m, d * d).sum(0)
    Qh, Qc, W = engine_quantities(p, body_after, Eh, Ec, d)
    eta = W / Qh
    J = a * s[0] - p[v[0]] * s[1]
    return dict(eta_star=es, eta_C=etaC, eta=eta, W=W, k=k, l=l, R=R,
                cat_err=float(np.abs(cat_after - cat_before).max()),
                formula_err=abs(eta - (1 - l * y / (k * x))), J=J)


def main():
    fails = []; n_ok = 0
    for d in (2, 3):
        for _ in range(12 if d == 3 else 40):
            t = trial(d)
            if t is None: continue
            if isinstance(t, tuple): fails.append(t); continue
            ok = (t["eta_star"] < t["eta_C"] - 1e-12 and t["cat_err"] < 1e-14 and t["W"] > 0 and t["J"] > 0
                  and t["eta"] > t["eta_star"] + 1e-12 and t["formula_err"] < 1e-9)
            n_ok += ok
            if not ok: fails.append(("trial", d, t))
            if d == 3 and n_ok % 4 == 1:
                print(f"d={d}: eta*={t['eta_star']:.5f} eta={t['eta']:.5f} eta_C={t['eta_C']:.5f} k={t['k']} l={t['l']} "
                      f"W={t['W']:.3e} catalyst err={t['cat_err']:.1e}")
    print(f"trials passing all checks: {n_ok}")
    # random 2-body unitaries never beat eta_* (d = 2)
    beat = 0
    for _ in range(200):
        d = 2; Eh = np.array([0, rng.uniform(0.5, 2)]); Ec = np.array([0, rng.uniform(0.5, 2)])
        bh = rng.uniform(0.1, 1); bc = bh * rng.uniform(1.5, 5); p = np.kron(gibbs(Eh, bh), gibbs(Ec, bc))
        es = eta_star(p, Eh, Ec, d)
        if es is None: continue
        Z = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4)); U, _ = np.linalg.qr(Z)
        q = (np.abs(U) ** 2) @ p
        Qh, Qc, W = engine_quantities(p, q, Eh, Ec, d)
        if W > 1e-12 and W / Qh > es + 1e-9: beat += 1
    print(f"random unitaries beating eta_*: {beat} (expected 0)")
    if beat: fails.append("unitary-beats-eta-star")
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:3]}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
