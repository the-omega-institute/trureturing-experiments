"""Numerical check of the counterexample family to Zietak's Conjecture 1.1 (Banach Center Publ. 112 (2017), DOI 10.4064/bc112-0-16).

Family (3x3): A = I, M = span{H, E33}, H = [[1,0,c],[0,-q,0],[c,0,0]], q > 1, c > 0.
Every residual A - B (B in M) is R(s,t) = [[1+s,0,cs],[0,1-qs,0],[cs,0,t]].
Strict spectral approximant: R(0,0) = diag(1,1,0), i.e. t = 0.
Claimed Schatten-p limit: t_p -> -1, i.e. B_p -> 2 E33 != E33.

Family (4x4): A = diag(1,1,3d/4,d/4), M = span{H (+) 0, diag(0,0,1,-1)}, 0 < d < 1.
Residual R(s,t) (+) [d - t]; strict residual t = d/2; claimed limit t_p -> d - 1.

The c_p objective is jointly convex in (s,t), so the minimum over t is convex in s and nested
golden-section search finds the real minimizer. Adjoint invariance makes the unique minimizer real;
a complex perturbation scan checks that no nearby complex point is lower.
"""
import numpy as np

bad = 0

def R3(s, t, q, c):
    return np.array([[1 + s, 0, c * s], [0, 1 - q * s, 0], [c * s, 0, t]], dtype=complex)

def R4(s, t, q, c, d):
    m = np.zeros((4, 4), dtype=complex)
    m[:3, :3] = R3(s, t, q, c)
    m[3, 3] = d - t
    return m

def logF(mat, p):
    sv = np.linalg.svd(mat, compute_uv=False)
    m = sv.max()
    return np.log(m) + np.log(np.sum((sv / m) ** p)) / p

def golden(f, a, b, it=160):
    g = (np.sqrt(5) - 1) / 2
    x1, x2 = b - g * (b - a), a + g * (b - a)
    f1, f2 = f(x1), f(x2)
    for _ in range(it):
        if f1 < f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - g * (b - a); f1 = f(x1)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + g * (b - a); f2 = f(x2)
    return (a + b) / 2

def minimize(build, p):
    def inner(s):
        t = golden(lambda t: logF(build(s, t), p), -2.5, 2.5)
        return t, logF(build(s, t), p)
    s = golden(lambda s: inner(s)[1], -1.0, 1.0)
    t, f = inner(s)
    return s, t, f

def complex_scan(build, s, t, f, p, rng):
    worst = np.inf
    for _ in range(200):
        ds = (rng.normal() + 1j * rng.normal()) * 1e-3
        dt = (rng.normal() + 1j * rng.normal()) * 1e-3
        worst = min(worst, logF(build(s + ds, t + dt), p) - f)
    return worst

rng = np.random.default_rng(15)
P = [10, 50, 200, 1000, 3000]

# Strict approximant face (3x3): ||R(s,t)||_inf <= 1 forces s = 0 (sampled check of the two disk constraints).
for q, c in [(2, 1), (1.5, 0.5), (3, 2), (1.1, 0.1)]:
    viol = 0
    for _ in range(20000):
        s = complex(rng.normal() * 0.3, rng.normal() * 0.3)
        if abs(s) < 1e-9:
            continue
        if abs(1 + s) <= 1 and abs(1 - q * s) <= 1:
            viol += 1
    if viol:
        print("FACE-VIOLATION", q, c, viol); bad += 1

print("3x3 family: strict residual diag(1,1,0) (t=0); expected t_p -> -1")
for q, c in [(2, 1), (1.5, 0.5), (3, 2), (1.1, 0.1)]:
    build = lambda s, t, q=q, c=c: R3(s, t, q, c)
    ts = []
    for p in P:
        s, t, f = minimize(build, p)
        gap = complex_scan(build, s, t, f, p, rng)
        ts.append(t)
        print(f"  q={q} c={c} p={p:5d} s={s:.3e} t={t:+.6f} ||R||_p={np.exp(f):.7f} complex_gap={gap:.2e}")
        if not (s > 0): print("  FAIL s_p > 0"); bad += 1
        if gap < -1e-12: print("  FAIL complex point lower"); bad += 1
        if not (np.exp(f) < 2 ** (1 / p)): print("  FAIL below strict residual"); bad += 1
    if not all(ts[i + 1] < ts[i] for i in range(len(ts) - 1)): print("  FAIL t_p monotone"); bad += 1
    if not (ts[-1] < -0.99): print("  FAIL t_3000 < -0.99"); bad += 1

print("4x4 family: strict residual t = d/2; expected t_p -> d - 1")
for q, c, d in [(2, 1, 0.5), (3, 2, 0.25)]:
    build = lambda s, t, q=q, c=c, d=d: R4(s, t, q, c, d)
    ts = []
    for p in P:
        s, t, f = minimize(build, p)
        ts.append(t)
        print(f"  q={q} c={c} d={d} p={p:5d} s={s:.3e} t={t:+.6f} target={d-1:+.3f}")
    if not (abs(ts[-1] - (d - 1)) < 0.02): print("  FAIL t_3000 near d-1"); bad += 1

print("bad=", bad)
raise SystemExit(1 if bad else 0)
