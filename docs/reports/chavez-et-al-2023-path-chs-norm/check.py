"""Chávez–Fullerton–LaFortune–Linarez–Liyanage–Son–Ting arXiv:2311.04689v3, Problem 1: compute ||P_n||_d, d >= 4.
Claimed closed form: ||P_n||_{2m}^{2m} = sum_{k=0}^{floor(m/(n+1))} [C(L, m-k(n+1)) - 2 C(L, m-k(n+1)-1) + C(L, m-k(n+1)-2)],
L = 2m+n-1, C(L,r) = 0 outside 0 <= r <= L.
Checks: (1) coefficients of 1/Q_n(x), Q_n = det(I - sqrt(x) A_n) (Q_0=Q_1=1, Q_n = Q_{n-1} - x Q_{n-2}), exact, n<=15, m<=25;
(2) h_{2m} of the actual path eigenvalues (numpy), n<=7, m<=6; (3) sample values; (4) control: dropping the
k>=1 wrap-around terms fails once m >= n+1.
"""
import sys
import numpy as np
import sympy as sp

x = sp.symbols('x')


def Q(n):
    a, b = sp.Integer(1), sp.Integer(1)
    if n <= 1: return sp.Integer(1)
    for _ in range(n - 1): a, b = b, sp.expand(b - x * a)
    return b


def C(L, r): return sp.binomial(L, r) if 0 <= r <= L else 0


def formula(n, m, kmax=None):
    L = 2 * m + n - 1
    K = m // (n + 1) if kmax is None else kmax
    return sum(C(L, m - k * (n + 1)) - 2 * C(L, m - k * (n + 1) - 1) + C(L, m - k * (n + 1) - 2) for k in range(K + 1))


def hcomplete(lams, d):
    c = np.zeros(d + 1); c[0] = 1
    for l in lams:
        new = np.zeros(d + 1)
        for i in range(d + 1): new[i] = sum(c[j] * l ** (i - j) for j in range(i + 1))
        c = new
    return c[d]


def main():
    fails = []; cnt = 0
    for n in range(1, 16):
        s = sp.series(1 / Q(n), x, 0, 26).removeO()
        for m in range(1, 26):
            cnt += 1
            if s.coeff(x, m) != formula(n, m): fails.append(("series", n, m))
    print(f"series coefficients compared: {cnt}, mismatches: {len(fails)}")
    for n in range(1, 8):
        A = np.diag(np.ones(n - 1), 1) + np.diag(np.ones(n - 1), -1)
        lam = np.linalg.eigvalsh(A)
        for m in range(1, 7):
            if abs(hcomplete(lam, 2 * m) - float(formula(n, m))) > 1e-6: fails.append(("eig", n, m))
    print("eigenvalue checks (n<=7, m<=6) done")
    samples = {(4, 2): 8, (4, 3): 21, (4, 4): 55, (3, 4): 16}
    for (n, m), v in samples.items():
        if formula(n, m) != v: fails.append(("sample", n, m))
    print("samples:", {k: int(formula(*k)) for k in samples})
    ctrl = sum(1 for n in range(1, 8) for m in range(n + 1, 20) if formula(n, m, kmax=0) != formula(n, m))
    print(f"control (no wrap-around terms) differs in {ctrl} cases (expected > 0)")
    if ctrl == 0: fails.append("control")
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:5]}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
