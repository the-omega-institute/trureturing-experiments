"""Petr-Turek, arXiv:2501.07269v2, Question 5.5 (are the eigenvalues lambda_1..lambda_k of
Theorem 1.5 all distinct?).

1. Evaluates the closed formulas of Theorem 1.5 (gcd(n,k) = 1) in exact integer arithmetic.
2. Validates them against the brute-force spectrum, with multiplicities, of the (n,k)-wreath
   matrix M(F,F') = |F cap F'| for (n,k) = (5,2), (7,2), (7,3), (8,3).
3. Reports the collision lambda_3 = lambda_5 at (n,k) = (14,5), the census of coprime collisions
   for n < 120, and checks the factorisation
   lambda_3(n,5) - lambda_5(n,5) = 60 (n-10)! (n-2)(n-7)(n-9)(n-10)(n-14) for 11 <= n < 200.
Requires numpy.  Usage: python3 check.py  (prints RESULT PASS / RESULT FAIL).
"""
from math import factorial as f, comb, gcd
from collections import Counter
import itertools
import numpy as np


def lam(n, k, l):
    if l == 1:
        g = gcd(n, k)
        if n % k:
            num, den = n * f(n - k) * f(k), 2 * g * f(g) ** (n // g)
        else:
            num, den = n * n * f(n - k), k * k * f(k) ** (n // k - 1) * f(n // k)
        assert num % den == 0
        return num // den
    s = sum((-1) ** (l - j) * f(k - j) * f(n - k - l + j) for j in range(l + 1))
    num = ((-1) ** l * f(l) * f(k) * f(n - k - l) * (n - 2 * k - 1) * comb(k, l)
           - f(k) * f(n - k) + 2 * f(l) * comb(k + 1, l + 1) * s)
    assert num % 2 == 0
    return num // 2


def wreaths(n, k):
    W = set()
    for p in itertools.permutations(range(1, n)):
        pi = (0,) + p
        W.add(frozenset(frozenset(pi[(i * k + t) % n] for t in range(k)) for i in range(n)))
    return list(W)


def main():
    ok = True
    for n, k in [(5, 2), (7, 2), (7, 3), (8, 3)]:
        W = wreaths(n, k)
        M = np.array([[len(a & b) for b in W] for a in W], dtype=float)
        spec = Counter(np.round(np.linalg.eigvalsh(M)).astype(int).tolist())
        pred = Counter({lam(n, k, 1): 1})
        for l in range(2, k + 1):
            pred[lam(n, k, l)] += comb(n, l) - comb(n, l - 1)
        pred[0] += len(W) - sum(pred.values())
        match = spec == pred
        ok &= match
        print(f"(n,k)=({n},{k}) wreaths={len(W)} spectrum={sorted(spec.items())} formula_match={match}")
    L = [lam(14, 5, l) for l in range(1, 6)]
    print("(n,k)=(14,5) lambda_1..lambda_5 =", L, " multiplicities", [1] + [comb(14, l) - comb(14, l - 1) for l in range(2, 6)])
    ok &= L[2] == L[4] == 17107200
    hits = [(n, k, [(i + 1, j + 1) for i in range(k) for j in range(i + 1, k) if lam(n, k, i + 1) == lam(n, k, j + 1)])
            for n in range(5, 120) for k in range(1, n) if 2 * k < n and gcd(n, k) == 1
            if len({lam(n, k, l) for l in range(1, k + 1)}) < k]
    print("coprime collisions for n < 120:", hits)
    fac = all(lam(n, 5, 3) - lam(n, 5, 5) == 60 * f(n - 10) * (n - 2) * (n - 7) * (n - 9) * (n - 10) * (n - 14)
              for n in range(11, 200) if gcd(n, 5) == 1)
    print("factorisation lambda_3(n,5) - lambda_5(n,5) holds for 11 <= n < 200, gcd(n,5) = 1:", fac)
    ok &= fac
    print("RESULT", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
