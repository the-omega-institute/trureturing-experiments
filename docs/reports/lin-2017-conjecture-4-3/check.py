"""Lin, Commun. Contemp. Math. 19 (2017) 1650044, Conjecture 4.3: exact refutation at p = 4.

Claim: for positive semidefinite A, B and every p > 0, lambda(A^2 + |BA|^p) majorizes
lambda(A^2 + |AB|^p), where |X| = (X* X)^(1/2).

For real symmetric A, B: |BA|^4 = (A B^2 A)^2 and |AB|^4 = (B A^2 B)^2.
"""
import sys
import numpy as np
import sympy as sp

R_ = sp.Rational
A = sp.diag(3, 1, R_(1, 2))
B = sp.Matrix([[R_(1, 2), R_(1, 4), R_(3, 4)], [R_(1, 4), R_(3, 4), 0], [R_(3, 4), 0, R_(5, 2)]])
ok = True

def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name)
    ok = ok and bool(cond)

lead = lambda M: [M[:k, :k].det() for k in range(1, M.shape[0] + 1)]
check("A positive definite (leading minors %s)" % lead(A), all(m > 0 for m in lead(A)))
check("B positive definite (leading minors %s)" % lead(B), all(m > 0 for m in lead(B)))

D = A * B * B * A          # |BA|^2
E = B * A * A * B          # |AB|^2
L = A * A + D * D          # A^2 + |BA|^4
R = A * A + E * E          # A^2 + |AB|^4
check("L, R symmetric", L == L.T and R == R.T)
check("trace L = trace R = %s" % L.trace(), L.trace() == R.trace())

theta = R_(7, 5)
mL = lead(L - theta * sp.eye(3))
check("L - 7/5 I positive definite (leading minors %s), so lambda_min(L) > 7/5" % mL, all(m > 0 for m in mL))
dR = (R - theta * sp.eye(3)).det()
check("det(R - 7/5 I) = %s < 0" % dR, dR < 0)
mR = lead(R)
check("R positive definite (leading minors %s)" % mR, all(m > 0 for m in mR))
# det(R - theta I) = prod (r_i - theta) < 0 forces an odd number of r_i below theta,
# so lambda_min(R) < 7/5 < lambda_min(L); with equal traces the sum of the two largest
# eigenvalues of R exceeds that of L.

lf = np.array(L.evalf(30).tolist(), dtype=float); rf = np.array(R.evalf(30).tolist(), dtype=float)
l = np.sort(np.linalg.eigvalsh(lf))[::-1]; r = np.sort(np.linalg.eigvalsh(rf))[::-1]
print("lambda(L) =", l, "partial sums", np.cumsum(l))
print("lambda(R) =", r, "partial sums", np.cumsum(r))
check("sum of two largest: L %.6f < R %.6f" % (l[:2].sum(), r[:2].sum()), l[:2].sum() < r[:2].sum())

# Consistency with Lin-Wolkowicz (2012) Corollary 2.3 at k = 2, which compares
# A^2 + (AB)^2 (BA)^2 with A^2 + (BA)^2 (AB)^2 (not |BA|^4 with |AB|^4).
LW1 = A * A + (A * B) ** 2 * (B * A) ** 2
LW2 = A * A + (B * A) ** 2 * (A * B) ** 2
a1 = np.sort(np.linalg.eigvals(np.array(LW1.evalf(30).tolist(), dtype=float)).real)[::-1]
a2 = np.sort(np.linalg.eigvals(np.array(LW2.evalf(30).tolist(), dtype=float)).real)[::-1]
check("Lin-Wolkowicz Cor. 2.3 (k = 2) holds here: partial sums %s >= %s" % (np.cumsum(a1), np.cumsum(a2)),
      all(np.cumsum(a1)[:2] >= np.cumsum(a2)[:2] - 1e-9) and abs(a1.sum() - a2.sum()) < 1e-9)
check("(AB)^2 (BA)^2 differs from |BA|^4", (A * B) ** 2 * (B * A) ** 2 != D * D)

# Exponent scan for this pair (floating point).
def pw(M, e):
    w, V = np.linalg.eigh(M); w = np.clip(w, 0, None); return (V * w ** e) @ V.T
af = np.array(A.evalf(30).tolist(), dtype=float); bf = np.array(B.evalf(30).tolist(), dtype=float)
df = af @ bf @ bf @ af; ef = bf @ af @ af @ bf
bad = []
for p in np.round(np.arange(0.25, 12.01, 0.25), 2):
    lp = np.sort(np.linalg.eigvalsh(af @ af + pw(df, p / 2)))[::-1]
    rp = np.sort(np.linalg.eigvalsh(af @ af + pw(ef, p / 2)))[::-1]
    g = (np.cumsum(lp) - np.cumsum(rp))[:-1].min() / lp[0]
    if g < -1e-12:
        bad.append(float(p))
print("exponents in 0.25..12 (step 0.25) where this pair violates the majorization:", bad)
check("this pair satisfies the majorization at p = 2", 2.0 not in bad)
print("RESULT", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
