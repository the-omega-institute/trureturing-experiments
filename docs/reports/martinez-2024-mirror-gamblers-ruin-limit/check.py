"""Martinez, arXiv:2412.07667v1, Conjecture 5.1 (mirror-step gambler's ruin, left (1-p)/2, right p, mirror (1-p)/2).

1. Solves the first-step recurrence f(x) = q f(x-1) + p f(x+1) + q f(N-x), f(0) = 0, f(N) = 1, exactly over Q
   for rational p and several N, and compares with the closed form
   F_N(x) = (1 - r^x)(1 - R r^(N-x)) / ((1 - R)(1 - r^N)),
   T = s/(ab) - 2, r = (T - sqrt(T^2 - 4))/2, K1 = 1 - a/r - b r, R = K1/c (a = c = q, b = p, s = a + b).
2. Checks K1 K2 = c^2 and R < 0, and that (1 - r)/(1 - R) equals the conjectured limit.
Requires sympy.  Usage: python3 check.py  (prints RESULT PASS / RESULT FAIL).
"""
import sympy as sp

def solve(N, p):
    q = (1 - p) / 2
    A = sp.zeros(N - 1, N - 1); b = sp.zeros(N - 1, 1)
    for x in range(1, N):
        i = x - 1; A[i, i] += 1
        for y, w in ((x - 1, q), (x + 1, p), (N - x, q)):
            if y == 0: continue
            if y == N: b[i] += w; continue
            A[i, y - 1] -= w
    return A.LUsolve(b)

def main():
    ok = True
    P = sp.symbols('p', positive=True)
    conj = (sp.sqrt((P + 1) * (1 - 3 * P + 4 * P ** 2)) - (1 - 2 * P) * (P + 1)) / (2 * P * (P + 1))
    a = (1 - P) / 2; b = P; c = (1 - P) / 2; s = a + b; T = s / (a * b) - 2
    r = (T - sp.sqrt(T ** 2 - 4)) / 2
    K1 = 1 - a / r - b * r; K2 = 1 - a * r - b / r; R = K1 / c
    for pv in [sp.Rational(1, 10), sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(7, 10), sp.Rational(9, 10)]:
        rr, RR = r.subs(P, pv), R.subs(P, pv)
        kk = sp.N((K1 * K2 - c ** 2).subs(P, pv), 40)
        lim = sp.N(((1 - r) / (1 - R)).subs(P, pv), 40); cj = sp.N(conj.subs(P, pv), 40)
        worst = 0
        for N in (3, 7, 12, 25):
            f = solve(N, pv)
            for x in range(1, N):
                F = (1 - rr ** x) * (1 - RR * rr ** (N - x)) / ((1 - RR) * (1 - rr ** N))
                worst = max(worst, abs(sp.N(f[x - 1] - F, 40)))
        f80 = sp.N(solve(80, pv)[0], 30)
        good = abs(kk) < 1e-30 and sp.N(RR) < 0 and 0 < sp.N(rr) < 1 and abs(lim - cj) < 1e-30 and worst < 1e-30
        ok &= bool(good)
        print(f"p={pv}: K1K2-c^2={float(kk):.1e} R={float(sp.N(RR)):.6f} r={float(sp.N(rr)):.6f} "
              f"max|f-F| (N=3,7,12,25)={float(worst):.1e} f_80(1)={float(f80):.15f} limit={float(cj):.15f} ok={good}")
    print("RESULT", "PASS" if ok else "FAIL")

if __name__ == "__main__":
    main()
