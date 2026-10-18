"""Büyükköroğlu–Aksoy–Dzhafarov, Commun. Optim. Theory 2024, Article 30, Conjecture 3.12: every robustly stable
segment [A1, A2] of real 2x2 matrices admits an affine parameter-dependent Lyapunov matrix (APDLM)
P(a) = (1-a) P1 + a P2 with P(a) > 0 and A(a)^T P(a) + P(a) A(a) < 0 for all a in [0, 1].
Numerical evidence (not a proof):
1. Random robustly stable segments (exact 2x2 test: trace < 0 at both ends, min_{[0,1]} det A(a) > 0) including
   near-boundary ones (min det close to 0): an LMI on a 401-point grid with margin finds P1, P2, which are then
   verified on a 5001-point grid. Count successes.
2. Control: segments that are not robustly stable are LMI-infeasible on the grid.
Uses cvxpy with the CLARABEL solver (venv /tmp/op-venv).
"""
import sys
import numpy as np
import cvxpy as cp

rng = np.random.default_rng(20261017)

def robust(A1, A2):
    if np.trace(A1) >= 0 or np.trace(A2) >= 0: return False, None
    a = np.linspace(0, 1, 2001)
    dets = [np.linalg.det((1 - x) * A1 + x * A2) for x in a]
    # exact min of the quadratic det on [0,1]
    D = lambda x: np.linalg.det((1 - x) * A1 + x * A2)
    c0 = D(0); c2 = (D(1) + D(-1)) / 2 - c0; c1 = (D(1) - D(-1)) / 2
    cands = [0.0, 1.0] + ([-c1 / (2 * c2)] if abs(c2) > 1e-14 and 0 < -c1 / (2 * c2) < 1 else [])
    m = min(D(x) for x in cands)
    return m > 0, m

def apdlm(A1, A2, grid=401, eps=1e-6):
    P1 = cp.Variable((2, 2), symmetric=True); P2 = cp.Variable((2, 2), symmetric=True); s = cp.Variable()
    cons = []
    for x in np.linspace(0, 1, grid):
        A = (1 - x) * A1 + x * A2; P = (1 - x) * P1 + x * P2
        cons += [P >> s * np.eye(2), A.T @ P + P @ A << -s * np.eye(2)]
    cons += [cp.trace(P1) + cp.trace(P2) <= 4]
    prob = cp.Problem(cp.Maximize(s), cons)
    try:
        prob.solve(solver='CLARABEL')
    except Exception:
        prob.solve(solver='SCS')
    if P1.value is None: return False, None
    P1v, P2v = P1.value, P2.value
    worst = np.inf
    for x in np.linspace(0, 1, 5001):
        A = (1 - x) * A1 + x * A2; P = (1 - x) * P1v + x * P2v
        worst = min(worst, np.linalg.eigvalsh(P).min(), -np.linalg.eigvalsh(A.T @ P + P @ A).max())
    return worst > 0, worst

def main():
    fails = []; ok = 0; tot = 0; hard = 0; minmargin = np.inf
    while tot < 300:
        A1 = rng.normal(size=(2, 2)) * rng.choice([1, 3]); A2 = rng.normal(size=(2, 2)) * rng.choice([1, 3])
        rs, m = robust(A1, A2)
        if not rs: continue
        tot += 1
        good, w = apdlm(A1, A2)
        ok += good; minmargin = min(minmargin, w if w is not None else -np.inf)
        if not good: fails.append(("random", A1.tolist(), A2.tolist(), m))
    print(f"random robustly stable segments: {ok}/{tot} admit a verified APDLM (smallest verified margin {minmargin:.2e})")
    # near-boundary segments: rescale the second endpoint towards losing robust stability
    nb_ok = 0; nb = 0
    while nb < 60:
        A1 = rng.normal(size=(2, 2)); A2 = rng.normal(size=(2, 2))
        rs, m = robust(A1, A2)
        if not rs or m > 0.05: continue
        nb += 1
        good, w = apdlm(A1, A2)
        nb_ok += good
        if not good: fails.append(("near-boundary", A1.tolist(), A2.tolist(), m))
    print(f"near-boundary segments (min det <= 0.05): {nb_ok}/{nb} admit a verified APDLM")
    # control: not robustly stable
    ctrl = 0; cn = 0
    while cn < 40:
        A1 = rng.normal(size=(2, 2)); A2 = rng.normal(size=(2, 2))
        if np.trace(A1) >= 0 or np.trace(A2) >= 0: continue
        rs, m = robust(A1, A2)
        if rs or np.linalg.eigvals(A1).real.max() >= 0 or np.linalg.eigvals(A2).real.max() >= 0: continue
        cn += 1
        good, w = apdlm(A1, A2)
        ctrl += (not good)
    print(f"control (both endpoints Hurwitz, segment not robustly stable): {ctrl}/{cn} infeasible (expected all)")
    if ctrl != cn: fails.append("control")
    print("RESULT", "PASS" if not fails else f"FAIL {len(fails)} {fails[:2]}")
    sys.exit(0 if not fails else 1)

if __name__ == "__main__":
    main()
