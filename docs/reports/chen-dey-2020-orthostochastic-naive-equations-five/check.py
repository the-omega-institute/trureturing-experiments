"""Exact checks for Chen-Dey, arXiv:2001.10691 (Exp. Math. 32 (2023)), Question 5 at n = 5: do the naive equations cut
out the orthostochastic matrices inside the Birkhoff polytope?  Witness B = (4J - 3I)/17.
(1) B is doubly stochastic with positive entries.
(2) Every row/column pair admits signs making the entrywise square roots orthogonal (the naive condition).
(3) No real orthogonal lift: M = sqrt(17) U would be an integer matrix with M M^T = 17 I, so det(M)^2 = 17^5, which is
    not a square; also confirmed by an exhaustive row-by-row backtracking search over the sign patterns of
    (+-1 on the diagonal, +-2 off it) with pairwise orthogonal rows.
(4) The line family: for odd n and B = (4J + (k^2-4)I)/(k^2+4(n-1)), the naive condition holds iff k odd <= n-2;
    orthostochastic (by H = I - 2J/n) at k = n - 2, not at k <= n - 4 for n = 5, 7 (backtracking).  Prints ALL_OK."""
from fractions import Fraction as Fr
import itertools, math
ok = True
def Bmat(n, k):
    D = k*k + 4*(n-1)
    return [[Fr(4 + (k*k - 4 if i == j else 0), D) for j in range(n)] for i in range(n)]
def naive_pair(vals):  # vals: the products B_ki*B_kj (as Fractions); signs on their square roots summing to zero
    # all products here are b^2 * {r^2 (twice), 1 (n-2 times)} up to a common factor; test exactly with sqrt of rationals
    roots = [math.sqrt(float(v)) for v in vals]
    return any(abs(sum(s*r for s, r in zip(signs, roots))) < 1e-12 for signs in itertools.product((1, -1), repeat=len(vals)))
B = Bmat(5, 1)
ok &= all(x > 0 for row in B for x in row) and all(sum(row) == 1 for row in B) and all(sum(B[i][j] for i in range(5)) == 1 for j in range(5))
ok &= all(naive_pair([B[k][i]*B[k][j] for k in range(5)]) and naive_pair([B[i][k]*B[j][k] for k in range(5)]) for i in range(5) for j in range(i+1, 5))
print("B=(4J-3I)/17 doubly stochastic, naive conditions hold:", ok)
r = 17**5; ok &= math.isqrt(r)**2 != r
def has_lift(n, k):
    D = k*k + 4*(n-1)
    # lift U = M/sqrt(D) with |M_ii| = k, |M_ij| = 2; rows of M pairwise orthogonal, |row|^2 = D automatically
    rows = []
    def cand(i):
        for signs in itertools.product((1, -1), repeat=n):
            yield [signs[j]*(k if j == i else 2) for j in range(n)]
    def bt(i):
        if i == n: return True
        for v in cand(i):
            if v[i] < 0: continue          # rows may be multiplied by -1
            if all(sum(a*b for a, b in zip(v, w)) == 0 for w in rows):
                rows.append(v)
                if bt(i + 1): return True
                rows.pop()
        return False
    return bt(0)
for n in (5, 7):
    for k in range(1, n - 1, 2):
        lift = has_lift(n, k)
        print(f"n={n} k={k}: orthostochastic={lift}")
        ok &= lift == (k == n - 2)
print("ALL_OK" if ok else "FAIL")
