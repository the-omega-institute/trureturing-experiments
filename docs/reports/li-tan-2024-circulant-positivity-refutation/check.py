"""Li-Tan, arXiv:2412.08022v1, Section III conjecture: H in H_n (PSD, unit diagonal) implies C_H PSD.

Witness H(a,e,f) = [[1,a,e,a],[a,1,a,f],[e,a,1,a],[a,f,a,1]] with (1+sqrt(ef))/2 < a < sqrt((1+e)(1+f))/2.
1. Exact rational member a = 5/8, e = 1/9, f = 4/9: H positive definite (Sylvester minors, exact), C_H computed
   from the source's z_k = (prod_l |h_{l,l+k}|)^(1/n) exp(i mean arg), first row (1, 5/8, 2/9, 5/8), and the
   vector (1,-1,1,-1) gives v* C_H v = -1/9 (eigenvalue -1/36).
2. The n = 3 case (proved by the authors) on random samples: C_H stays PSD.
3. Random n = 4 samples: frequency of indefinite C_H, and the family gap check.
Requires numpy, sympy.  Usage: python3 check.py  (prints RESULT PASS / RESULT FAIL).
"""
import numpy as np, sympy as sp

def CH(H):
    n = H.shape[0]
    z = []
    for k in range(n):
        ent = [H[l, (l + k) % n] for l in range(n)]
        mod = np.prod([abs(e) for e in ent]) ** (1 / n)
        ph = sum(np.angle(e) for e in ent) / n
        z.append(mod * np.exp(1j * ph))
    C = np.zeros((n, n), complex)
    for i in range(n):
        for j in range(n):
            C[i, j] = z[(j - i) % n]
    return C

def main():
    ok = True
    a, e, f = sp.Rational(5, 8), sp.Rational(1, 9), sp.Rational(4, 9)
    H = sp.Matrix([[1, a, e, a], [a, 1, a, f], [e, a, 1, a], [a, f, a, 1]])
    minors = [H[:k, :k].det() for k in range(1, 5)]
    print("leading minors of H:", minors)
    ok &= all(m > 0 for m in minors)
    z2 = sp.root(e * f * e * f, 4)
    print("z = (1, 5/8,", z2, ", 5/8)")
    C = sp.Matrix(4, 4, lambda i, j: [1, a, z2, a][(j - i) % 4])
    v = sp.Matrix([1, -1, 1, -1])
    val = sp.simplify((v.T * C * v)[0])
    print("v^T C_H v =", val, " eigenvalue", sp.simplify(val / 4))
    ok &= val < 0 and z2 == sp.Rational(2, 9)
    Cn = CH(np.array(H.tolist(), dtype=float))
    print("numeric C_H min eigenvalue", np.linalg.eigvalsh(Cn).min())
    rng = np.random.default_rng(0)
    bad3 = bad4 = 0
    for _ in range(20000):
        for n in (3, 4):
            G = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            V = G / np.linalg.norm(G, axis=0)
            Hs = V.conj().T @ V
            m = np.linalg.eigvalsh((CH(Hs) + CH(Hs).conj().T) / 2).min()
            if m < -1e-10:
                if n == 3: bad3 += 1
                else: bad4 += 1
    print(f"random Gram samples (20000 each): indefinite C_H for n=3: {bad3}, n=4: {bad4}")
    ok &= bad3 == 0
    print("RESULT", "PASS" if ok else "FAIL")

if __name__ == "__main__":
    main()
