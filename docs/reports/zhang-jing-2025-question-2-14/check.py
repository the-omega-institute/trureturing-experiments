"""Zhang-Jing, ELA 41 (2025) 595-601, Question 2.14: is [[f(|A|), A*], [A, g(|A|)]] >= 0 iff A is normal, for f, g
with f g = x^2 (x >= 0), f >= g (x >= 1), g >= f (0 <= x <= 1)?

Counterexample family: g(x) = 1 + (x-1)/(1+(x-1)^2), f(x) = x^2/g(x), A_s = [[0, 1+1/s], [1+s, 0]] (s > 1).
Exact check at s = 2 (Fractions): |A| = diag(3, 3/2), g(|A|) = (7/5) I, c*M = Q*Q with Q = [A, cI], A nonnormal.
Numerical checks (numpy): admissibility of f, g on a grid; M_s >= 0 and A_s nonnormal for several s; for the injective
pair f = x^(3/2), g = x^(1/2) (a = 3/4 in the open range 1/2 <= a <= 1 of the source's Theorem 2.3) the same A_s gives M not PSD.
Usage: python3 check.py
"""
import sys
from fractions import Fraction as Fr
import numpy as np

def g(x): return 1 + (x - 1) / (1 + (x - 1) ** 2)
def f(x): return x * x / g(x)

def mat(a): return [[Fr(x) for x in r] for r in a]
def mul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def T(A): return [list(r) for r in zip(*A)]

def main():
    ok = True
    # admissibility on a grid
    xs = np.linspace(0, 20, 20001)
    adm = all(abs(f(x) * g(x) - x * x) < 1e-9 for x in xs) and all(f(x) >= g(x) - 1e-12 for x in xs if x >= 1) and all(g(x) >= f(x) - 1e-12 for x in xs if x <= 1) and min(g(x) for x in xs) > 0
    print(f"f g = x^2, f >= g on [1, 20], g >= f on [0, 1], g > 0 (grid): {adm}"); ok &= adm
    # exact s = 2
    s = Fr(2); a, b = 1 + s, 1 + 1 / s; c = g(a)
    A = mat([[0, b], [a, 0]]); D = mat([[a, 0], [0, b]])
    AtA = mul(T(A), A); AAt = mul(A, T(A))
    exact = AtA == mul(D, D) and AAt != AtA and g(a) == g(b) == Fr(7, 5)
    M = [[a * a / c, 0, 0, a], [0, b * b / c, b, 0], [0, b, c, 0], [a, 0, 0, c]]  # [[f(D), A*], [A, g(D)]], A* = A^T
    Q = [[0, b, c, 0], [a, 0, 0, c]]
    cert = mul(T(Q), Q) == [[c * x for x in r] for r in M]
    print(f"s = 2: |A| = diag(3, 3/2) {AtA == mul(D, D)}, nonnormal {AAt != AtA}, g(3) = g(3/2) = 7/5 {g(a) == g(b) == Fr(7, 5)}; c M = Q^T Q {cert}"); ok &= exact and cert
    # numerical family
    fam = True
    for sv in [1.5, 2, 3, 5, 10, 100]:
        a, b = 1 + sv, 1 + 1 / sv; Dn = np.diag([a, b]); An = np.array([[0, b], [a, 0]])
        Mn = np.block([[np.diag([f(a), f(b)]), An.T], [An, np.diag([g(a), g(b)])]])
        fam &= np.linalg.eigvalsh(Mn).min() > -1e-9 and not np.allclose(An.T @ An, An @ An.T)
    print(f"family s in {{1.5, 2, 3, 5, 10, 100}}: M >= 0 and A nonnormal: {fam}"); ok &= fam
    # contrast: injective pair
    a, b = 3.0, 1.5; An = np.array([[0, b], [a, 0]])
    Mp = np.block([[np.diag([a ** 1.5, b ** 1.5]), An.T], [An, np.diag([a ** 0.5, b ** 0.5])]])
    mn = np.linalg.eigvalsh(Mp).min()
    print(f"contrast f = x^(3/2), g = x^(1/2) at s = 2: min eigenvalue {mn:.6f} (< 0)"); ok &= mn < 0
    print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)

main()
