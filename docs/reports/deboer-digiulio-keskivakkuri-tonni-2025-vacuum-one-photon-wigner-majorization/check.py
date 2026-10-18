"""Wigner majorization along the vacuum/one-photon family of de Boer, Di Giulio,
Keski-Vakkuri and Tonni, arXiv:2412.19698 (PRA 112, 032405, 2025), Section V.A.1.

State family, Eqs. (100)-(101):  W_u(x,p) = (1-u) W_0 + u W_1
  W_n(x,p) = (-1)^n / pi * exp(-(x^2+p^2)) * L_n(2(x^2+p^2)),  so
  W_u(x,p) = exp(-s) (1 - 2u + 2u s) / pi,  s = x^2 + p^2.
Proposal 1 / Eqs. (98)-(99):  rho_v >=_W rho_u  iff  I_t(v) >= I_t(u) for all t >= 0,
  I_t(u) = integral over R^2 of max(|W_u| - t, 0).

Claim under test (preregistered):
  (A) 1/2 <= u <= v <= 1  ==>  I_t(u) <= I_t(v)   for all t >= 0;
  (B) 0 <= u <= v <= 1/2  ==>  I_t(v) <= I_t(u)   for all t >= 0.
Auxiliary closed forms used by the proof route:
  I_0(u) = 1 for u <= 1/2,  I_0(u) = 4u exp(-1 + 1/(2u)) - 1 for u > 1/2,
  M(u) = int_0^inf f_u(s)^2 ds = 1/2 - u + u^2   (f_u = pi W_u in the variable s),
and the cross-branch incomparability for 0 < u < 1/2 < v <= 1-u.
Negative controls: the full family [0,1] is not a chain; the reversed direction on
each branch fails; a 2-D brute-force quadrature fixes the polar factor pi.
"""
import itertools, math, sys
import mpmath as mp
import numpy as np

mp.mp.dps = 30


def f(u, s):
    return mp.e**(-s) * (1 - 2*u + 2*u*s)


def bisect(g, a, b, it=120):
    ga = g(a)
    for _ in range(it):
        m = (a + b) / 2
        gm = g(m)
        if (gm > 0) == (ga > 0):
            a, ga = m, gm
        else:
            b = m
    return (a + b) / 2


def crossings(u, T):
    """Points where |f_u(s)| = T.  With a = 1-2u, b = 2u: f = e^{-s}(a + b s) is monotone on
    [0, s*] and [s*, inf), s* = 1 - a/b; if a < 0 it also changes sign at s0 = -a/b."""
    a, b = 1 - 2*u, 2*u
    g = lambda s: abs(f(u, s)) - T
    knots = [mp.mpf(0)]
    if b > 0:
        if a < 0:
            knots.append(-a / b)
        sstar = 1 - a / b
        if sstar > 0:
            knots.append(sstar)
    knots.append(mp.mpf(200))
    pts = []
    for lo, hi in zip(knots, knots[1:]):
        if (g(lo) > 0) != (g(hi) > 0):
            pts.append(bisect(g, lo, hi))
    if a < 0:
        pts.append(-a / b)
    return sorted(set(pts)), g


def I(u, t):
    """I_t(u) = pi * int_0^inf max(|f_u(s)|/pi - t, 0) ds = int_0^inf max(|f_u| - pi t, 0) ds."""
    u, T = mp.mpf(u), mp.pi * mp.mpf(t)
    pts, g = crossings(u, T)
    knots = [mp.mpf(0)] + pts + [mp.inf]
    total = mp.mpf(0)
    for a, b in zip(knots, knots[1:]):
        mid = a + 1 if b == mp.inf else (a + b) / 2
        if g(mid) > 0:
            total += mp.quad(lambda s: abs(f(u, s)) - T, [a, b])
    return total


def M(u):
    u = mp.mpf(u)
    return mp.quad(lambda s: f(u, s)**2, [0, mp.inf])


def brute_2d(u, t, L=7.0, n=2801):
    xs = np.linspace(-L, L, n)
    X, P = np.meshgrid(xs, xs)
    S = X**2 + P**2
    W = np.exp(-S) * (1 - 2*u + 2*u*S) / math.pi
    h = xs[1] - xs[0]
    return float(np.sum(np.maximum(np.abs(W) - t, 0.0)) * h * h)


def main():
    fails = []
    tmax = 1 / math.pi                      # max |W_u| <= 1/pi, so I_t = 0 beyond
    ts = [tmax * k / 40 for k in range(0, 41)]
    upper = [0.5 + k / 40 for k in range(0, 21)]
    lower = [k / 40 for k in range(0, 21)]
    table = {}
    for u in sorted(set(upper + lower + [0.1, 0.3, 0.6, 0.7, 0.9])):
        table[u] = [I(u, t) for t in ts]
    tol = mp.mpf('1e-18')
    # (A) increasing chain on [1/2, 1]
    nA = 0
    for u, v in itertools.combinations(upper, 2):
        for k in range(len(ts)):
            nA += 1
            if table[v][k] - table[u][k] < -tol:
                fails.append(('A', u, v, ts[k], table[u][k], table[v][k]))
    # (B) decreasing chain on [0, 1/2]
    nB = 0
    for u, v in itertools.combinations(lower, 2):
        for k in range(len(ts)):
            nB += 1
            if table[u][k] - table[v][k] < -tol:
                fails.append(('B', u, v, ts[k], table[u][k], table[v][k]))
    print(f"(A) pairs*thresholds checked {nA}; (B) checked {nB}; violations {len(fails)}")
    # closed forms
    cf = 0
    for u in table:
        exp0 = 1 if u <= 0.5 else 4*u*mp.e**(-1 + 1/(2*mp.mpf(u))) - 1
        e1 = abs(table[u][0] - exp0)
        e2 = abs(M(u) - (mp.mpf(1)/2 - u + mp.mpf(u)**2))
        cf = max(cf, e1, e2)
    print(f"closed forms I_0 and M: max abs error {mp.nstr(cf, 5)}")
    if cf > mp.mpf('1e-15'):
        fails.append(('closed-form', cf))
    # cross-branch incomparability 0 < u < 1/2 < v <= 1-u
    inc = 0
    for u in [k / 40 for k in range(1, 20)]:
        for v in [0.5 + k / 40 for k in range(1, 21)]:
            if v > 1 - u + 1e-12:
                continue
            Iu = table.get(u) or [I(u, t) for t in ts]
            Iv = table.get(v) or [I(v, t) for t in ts]
            up = any(Iv[k] - Iu[k] > tol for k in range(len(ts)))
            dn = any(Iu[k] - Iv[k] > tol for k in range(len(ts)))
            inc += 1
            if not (up and dn):
                fails.append(('incomparable', u, v, up, dn))
    print(f"cross-branch pairs with 0<u<1/2<v<=1-u checked {inc}")
    # negative controls
    c1 = any(table[0.1][k] - table[0.6][k] > tol for k in range(len(ts))) and \
         any(table[0.6][k] - table[0.1][k] > tol for k in range(len(ts)))
    revA = all(table[0.6][k] >= table[0.9][k] - tol for k in range(len(ts)))   # reversed (A) must fail
    revB = all(table[0.3][k] >= table[0.1][k] - tol for k in range(len(ts)))   # reversed (B) must fail
    print(f"control full-family chain fails (0.1 vs 0.6 cross): {c1}")
    print(f"control reversed (A) for (0.6,0.9) holds: {revA} (expected False)")
    print(f"control reversed (B) for (0.1,0.3) holds: {revB} (expected False)")
    if not c1:
        fails.append(('control-full-family', c1))
    if revA:
        fails.append(('control-reversed-A', revA))
    if revB:
        fails.append(('control-reversed-B', revB))
    rA = all(table[0.9][k] >= table[0.6][k] - tol for k in range(len(ts)))
    rB = all(table[0.1][k] >= table[0.3][k] - tol for k in range(len(ts)))
    strict_ok = any(table[0.9][k] > table[0.6][k] + tol for k in range(len(ts))) and \
        any(table[0.1][k] > table[0.3][k] + tol for k in range(len(ts)))
    print(f"strictness witnesses present: {strict_ok}; chains (0.6<=0.9),(0.3<=0.1): {rA},{rB}")
    # 2-D brute force against the radial reduction (polar factor pi)
    worst = 0.0
    for u, t in [(0.3, 0.0), (0.3, 0.1), (0.7, 0.0), (0.7, 0.05), (1.0, 0.12)]:
        b, r = brute_2d(u, t), float(I(u, t))
        worst = max(worst, abs(b - r))
        print(f"2-D brute force u={u} t={t}: grid {b:.8f} radial {r:.8f}")
    if worst > 1e-5:
        fails.append(('brute-2d', worst))
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:5]}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
