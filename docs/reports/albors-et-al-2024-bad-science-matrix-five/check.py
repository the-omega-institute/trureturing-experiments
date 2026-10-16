"""Checks for Albors-Bhatti-Ganjoo-Guo-Kunisky-Mukherjee-Stepin-Zeng, arXiv:2408.00933 (Involve 19 (2026)), Section 1.4:
is the displayed 5x5 matrix optimal, i.e. beta(A) <= (2 + 3 sqrt 3)/4 for every real 5x5 A with unit rows,
beta(A) = 2^-5 sum_{x in {-1,1}^5} ||A x||_inf ?
(1) beta(A*) = (2 + 3 sqrt 3)/4 exactly (sympy).
(2) Lemmas 2-3: for every antipodal-free S of the 5-cube with |S| = m <= 7, ||sum S||^2 <= L(m)^2,
    L(m) = (2 sqrt3 - 2) m + 8 - 4 sqrt 3 (exhaustive over orientations of m-subsets of x_1 = 1); and the
    Hamming bound q <= m^2 + 4m + 4 floor(m^2/4) is attained for m = 2, 4.
(3) Lemma 1 sanity for m >= 8 (random antipodal-free subsets) and random/locally optimized A never exceed the bound.
Prints ALL_OK.  Needs numpy and sympy."""
import itertools, math, random
import sympy as sp
import numpy as np
ok = True
r3 = sp.sqrt(3)
A = sp.Matrix([[1, 1, 0, 0, 1], [-1, 1, 0, 1, 0], [-1, 0, 0, -1, 1]]) / r3
A = A.col_join(sp.Matrix([[0, -1, 1, 1, 1], [0, 1, 1, -1, -1]]) / 2)
ok &= all(sp.simplify(sum(A[i, j]**2 for j in range(5)) - 1) == 0 for i in range(5))
cube = list(itertools.product((-1, 1), repeat=5))
vals = [max(abs(sp.nsimplify(sum(A[i, j]*x[j] for j in range(5)))) for i in range(5)) for x in cube]
beta = sp.nsimplify(sum(vals) / 32)
print("beta(A*) =", sp.simplify(beta), " target (2+3 sqrt3)/4:", sp.simplify(beta - (2 + 3*r3)/4) == 0,
      " counts:", {str(sp.simplify(v)): vals.count(v) for v in set(vals)})
ok &= sp.simplify(beta - (2 + 3*r3)/4) == 0
L = lambda m: (2*math.sqrt(3) - 2)*m + 8 - 4*math.sqrt(3)
R = [x for x in cube if x[0] == 1]
for m in range(0, 8):
    best = 0
    for sub in itertools.combinations(R, m):
        for signs in itertools.product((1, -1), repeat=m):
            u = [sum(s*x[k] for s, x in zip(signs, sub)) for k in range(5)]
            best = max(best, sum(c*c for c in u))
    hb = m*m + 4*m + 4*(m*m//4)
    print(f"m={m}: max ||u||^2 = {best}, Hamming bound {hb}, L(m)^2 = {L(m)**2:.4f}")
    ok &= best <= L(m)**2 + 1e-9 and best <= hb
    if m in (2, 4): ok &= best == hb and abs(best - L(m)**2) < 1e-9
rng = random.Random(3)
for _ in range(20000):
    m = rng.randint(8, 16)
    sub = rng.sample(R, m)
    signs = [rng.choice((1, -1)) for _ in sub]
    u = [sum(s*x[k] for s, x in zip(signs, sub)) for k in range(5)]
    ok &= sum(c*c for c in u) <= 16*m and sum(c*c for c in u) <= L(m)**2 + 1e-9
X = np.array(cube, dtype=float)
bound = (2 + 3*math.sqrt(3))/4
def beta_np(M): return np.abs(X @ M.T).max(axis=1).mean()
nprng = np.random.default_rng(5); best = 0.0
for _ in range(400):
    M = nprng.normal(size=(5, 5)); M /= np.linalg.norm(M, axis=1, keepdims=True)
    b = beta_np(M); step = 0.3
    for it in range(400):   # crude random local ascent
        N = M + step*nprng.normal(size=(5, 5)); N /= np.linalg.norm(N, axis=1, keepdims=True)
        bn = beta_np(N)
        if bn > b: M, b = N, bn
        else: step *= 0.995
    best = max(best, b)
print(f"local-ascent best beta over 400 starts = {best:.9f}, bound = {bound:.9f}")
ok &= best <= bound + 1e-9
print("ALL_OK" if ok else "FAIL")
