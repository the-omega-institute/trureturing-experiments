"""Fokkink–Joshi arXiv:2507.16956v4, Section 4 (after eq. (14)): the generalisation of Conjecture 20 says that
floor(s_n) = floor(alpha n + beta) for the tails s_n = u n + v + n^2/s_{n+1} of the continued fraction (14),
alpha the u-th metallic mean, beta = (v alpha^2 - alpha)/(alpha^2 + 1). For u = 2, v = 1 (alpha = 1 + sqrt2, beta = 1/2) it
fails at the Pell indices n_r = 1, 6, 35, 204, 1189, ... (n_{r+2} = 6 n_{r+1} - n_r).
Tails are computed by backward recursion from s_N = u N + v with N = n + 400 at 80 significant digits (the map
y -> u n + v + n^2/y contracts by a factor < 1/4, so truncation error is far below the printed precision).
Controls: the source's Theorem 24 case u = 1, v = phi - 1 has no failure for n <= 3000; the original Conjecture 20
(v = -1) for u = 1, 2, 3 is also tested for n <= 3000 (for Triage).
"""
import sys
from mpmath import mp, mpf, sqrt, floor
mp.dps = 80

def tails(u, v, n_max, extra=400):
    N = n_max + extra
    s = u * N + v
    out = {}
    for n in range(N - 1, 0, -1):
        s = u * n + v + mpf(n) ** 2 / s
        if n <= n_max: out[n] = s
    return out

def check(u, v, n_max):
    a = (u + sqrt(u * u + 4)) / 2; b = (v * a * a - a) / (a * a + 1)
    s = tails(u, v, n_max)
    fails = [n for n in range(1, n_max + 1) if floor(s[n]) != floor(a * n + b)]
    err = max(abs(s[n] - a * n - b) * n for n in range(1, n_max + 1))
    return fails, a, b, err, s

def main():
    bad = []
    fails, a, b, err, s = check(2, 1, 1300)
    pell = [1, 6]
    while pell[-1] < 1300: pell.append(6 * pell[-1] - pell[-2])
    pell = [p for p in pell if p <= 1300]
    print(f"u=2, v=1: alpha = {mp.nstr(a, 15)}, beta = {mp.nstr(b, 15)}; failures for n <= 1300: {fails}; Pell indices {pell}")
    for n in pell[:4]:
        print(f"   n={n}: s_n = {mp.nstr(s[n], 20)}, floor(s_n) = {int(floor(s[n]))}, floor(alpha n + beta) = {int(floor(a*n+b))}, n(s_n - alpha n - beta) = {mp.nstr(n*(s[n]-a*n-b), 10)}")
    if not set(pell) <= set(fails): bad.append("pell")
    if any(not (1 / mpf(8 * n) < s[n] - a * n - b < a / (4 * n)) for n in range(1, 1301)): bad.append("bounds")
    phi = (1 + sqrt(5)) / 2
    f24, *_ = check(1, phi - 1, 3000)
    print(f"control (Theorem 24) u=1, v=phi-1: failures for n <= 3000: {len(f24)}")
    if f24: bad.append("control")
    for u in (1, 2, 3):
        fc, a2, b2, e2, _ = check(u, -1, 3000)
        print(f"Triage: Conjecture 20 (v = -1), u={u}: failures for n <= 3000: {len(fc)} {fc[:6]}")
    print("RESULT", "PASS" if not bad else f"FAIL {bad}")
    sys.exit(0 if not bad else 1)

if __name__ == "__main__":
    main()
