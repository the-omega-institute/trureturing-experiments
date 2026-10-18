"""Glen–Simpson–Smyth arXiv:1701.00928 (EJC 24(3) P3.28), Open Problem 2: asymptotics of ET and ED.

Exact (paper's (3) and (6)), with a_m = L_sigma(m) sigma^{-m}, L_sigma(m) = (1/m) sum_{d|m} mu(d) sigma^{m/d}:
  ET = sum_{m=1}^n (n-m+1) a_m
  ED = sum_{m=1}^n L_sigma(m) sum_{s=1}^{floor(n/m)} (-1)^{s+1} C(n - s m + s, s) sigma^{-s m}
Claimed (r354d):
  ET = n ln n + (gamma - 1 + A) n + ln n + gamma + 1/2 + A - B + O(1/n),
  ED = n ln n - n ln ln n + (ln ln sigma - 1) n + O(n / ln n),
  A = sum_m (a_m - 1/m), B = sum_m m (a_m - 1/m).
Checks: Table 2 and Table 3 of the paper reproduced; n * |ET - formula| bounded; |ED - formula| / (n/ln n) bounded.
"""
import math, sys
import mpmath as mp

mp.mp.dps = 700


def mobius(k):
    res, p, x = 1, 2, k
    while p * p <= x:
        if x % p == 0:
            x //= p
            if x % p == 0:
                return 0
            res = -res
        p += 1
    return -res if x > 1 else res


def L(sig, m):
    return sum(mobius(d) * sig ** (m // d) for d in range(1, m + 1) if m % d == 0) // m


def ET(sig, n):
    return mp.fsum((n - m + 1) * mp.mpf(L(sig, m)) / mp.mpf(sig) ** m for m in range(1, n + 1))


def ED(sig, n):
    tot = mp.mpf(0)
    for m in range(1, n + 1):
        inner = mp.fsum((-1) ** (s + 1) * mp.binomial(n - s * m + s, s) * mp.mpf(sig) ** (-s * m) for s in range(1, n // m + 1))
        tot += L(sig, m) * inner
    return tot


def consts(sig, M=600):
    a = [mp.mpf(L(sig, m)) / mp.mpf(sig) ** m for m in range(1, M + 1)]
    A = mp.fsum(a[m - 1] - mp.mpf(1) / m for m in range(1, M + 1))
    B = mp.fsum(m * (a[m - 1] - mp.mpf(1) / m) for m in range(1, M + 1))
    return A, B


def main():
    fails = []
    t2 = {(2, 10): 18.62, (5, 10): 21.13, (2, 5): 7.31, (5, 5): 8.24}
    t3 = {(2, 10): 8.83, (5, 10): 15.13, (10, 30): 70.11, (20, 30): 78.33, (2, 30): 37.70}
    for (s, n), v in t2.items():
        if abs(float(ET(s, n)) - v) > 0.006: fails.append(("table2", s, n, float(ET(s, n))))
    for (s, n), v in t3.items():
        if abs(float(ED(s, n)) - v) > 0.006: fails.append(("table3", s, n, float(ED(s, n))))
    print("tables 2 and 3 reproduced:", not fails)
    g = mp.euler
    for sig in (2, 3, 5):
        A, B = consts(sig)
        print(f"sigma={sig}: A={mp.nstr(A, 12)} B={mp.nstr(B, 12)}")
        for n in (100, 1000, 10000):
            et = ET(sig, n)
            f = n * mp.log(n) + (g - 1 + A) * n + mp.log(n) + g + mp.mpf(1) / 2 + A - B
            err = (et - f) * n
            print(f"  ET n={n}: n*(ET - formula) = {mp.nstr(err, 8)}")
            if abs(err) > 50: fails.append(("ET", sig, n, float(err)))
        for n in (50, 200, 800):
            ed = ED(sig, n)
            f = n * mp.log(n) - n * mp.log(mp.log(n)) + (mp.log(mp.log(sig)) - 1) * n
            ratio = (ed - f) / (n / mp.log(n))
            print(f"  ED n={n}: ED={mp.nstr(ed, 10)} (ED - formula)/(n/ln n) = {mp.nstr(ratio, 6)}")
            if abs(ratio) > 20: fails.append(("ED", sig, n, float(ratio)))
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:4]}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
