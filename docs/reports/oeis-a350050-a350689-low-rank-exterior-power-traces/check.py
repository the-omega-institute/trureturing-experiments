"""OEIS A350689 and A350050 (S. Spezia): traces of exterior powers of three structured matrices.

Conjectures (each "trace of the k-th exterior power" is, as the entries state, up to sign the
coefficient of x^(n-k) of the characteristic polynomial):
- A350689: for M[i,j] = floor((j - i + 1)/2), |e_3(M)| = a(n) = n(1 - (-1)^n - 2(3 + (-1)^n)n^2 + 2n^4)/384 and
  e_k(M) = 0 for k > 3.
- A350050: for A[i,j] = n - |(n+1)/2 - j| - |(n+1)/2 - i| (A349107) and for
  B[i,j] = (n mod 2) + |(n+1)/2 - j| + |(n+1)/2 - i| (A349108), |e_2| = a(n) = (2n^4 - 6(-1)^n n^2 - 2n^2 + 3(-1)^n - 3)/96
  and e_k = 0 for k > 2.

Proof checked here: M = U V^T with n x 3 factors (M[i,j] = (j - i)/2 + 1/4 - (-1)^(i+j)/4), and A, B = U V^T with
n x 2 factors (constants and g_i = |(n+1)/2 - i|). Then charpoly(U V^T) = x^(n-r) charpoly(V^T U); the script checks
the factorizations, the closed forms of det(V^T U) (and e_2 for A, B) against the formulas, and the conjectures
themselves via exact characteristic polynomials for n <= N.
"""
import sys
from sympy import Matrix, Rational, floor, symbols, Abs, eye, zeros
N = int(sys.argv[1]) if len(sys.argv) > 1 else 22
x = symbols('x')
def a689(n): return Rational(n * (1 - (-1) ** n - 2 * (3 + (-1) ** n) * n ** 2 + 2 * n ** 4), 384)
def a050(n): return Rational(2 * n ** 4 - 6 * (-1) ** n * n ** 2 - 2 * n ** 2 + 3 * (-1) ** n - 3, 96)
def coeffs(Mx):
    p = Mx.charpoly(x).all_coeffs()   # p[k] = (-1)^k e_k
    return [((-1) ** k) * p[k] for k in range(len(p))]
ok = True
for n in range(2, N + 1):
    I = range(1, n + 1)
    M = Matrix(n, n, lambda i, j: floor(Rational((j + 1) - (i + 1) + 1, 2)))
    U = Matrix(n, 3, lambda i, a: [1, i + 1, (-1) ** (i + 1)][a])
    V = Matrix(n, 3, lambda j, a: [Rational(j + 1, 2) + Rational(1, 4), Rational(-1, 2), Rational(-((-1) ** (j + 1)), 4)][a])
    f689 = (U * V.T == M)
    K = V.T * U
    e = coeffs(M)
    r689 = (abs(e[3]) == a689(n)) if n >= 3 else (a689(n) == 0)
    z689 = all(e[k] == 0 for k in range(4, n + 1))
    k689 = (n < 3) or (abs(K.det()) == a689(n))
    c = Rational(n + 1, 2)
    g = [abs(c - i) for i in I]
    A = Matrix(n, n, lambda i, j: n - g[j] - g[i])
    B = Matrix(n, n, lambda i, j: (n % 2) + g[j] + g[i])
    UA = Matrix(n, 2, lambda i, a: [1, g[i]][a]); VA = Matrix(n, 2, lambda j, a: [n - g[j], -1][a])
    UB = Matrix(n, 2, lambda i, a: [1, g[i]][a]); VB = Matrix(n, 2, lambda j, a: [(n % 2) + g[j], 1][a])
    fA = (UA * VA.T == A); fB = (UB * VB.T == B)
    eA = coeffs(A); eB = coeffs(B)
    rA = abs(eA[2]) == a050(n); zA = all(eA[k] == 0 for k in range(3, n + 1))
    rB = abs(eB[2]) == a050(n); zB = all(eB[k] == 0 for k in range(3, n + 1))
    kA = abs((VA.T * UA).det()) == a050(n); kB = abs((VB.T * UB).det()) == a050(n)
    row = [f689, r689, z689, k689, fA, rA, zA, kA, fB, rB, zB, kB]
    if not all(row): ok = False; print('FAIL', n, row)
print(f"n = 2..{N}: factorizations, |e_3| = A350689(n), e_k = 0 (k > 3), |det V^T U| = A350689(n); A349107/A349108 matrices: factorizations, |e_2| = A350050(n), e_k = 0 (k > 2), |det V^T U| = A350050(n): {ok}")
D689 = [0, 0, 0, 1, 4, 15, 36, 84, 160, 300, 500, 825, 1260, 1911, 2744, 3920, 5376, 7344]
D050 = [0, 0, 0, 2, 4, 14, 24, 52, 80, 140, 200, 310, 420, 602, 784, 1064, 1344, 1752, 2160]
ok_data = [a689(n) for n in range(0, len(D689))] == D689 and [a050(n) for n in range(0, len(D050))] == D050
print(f"published data (first terms): {ok_data}")
good = ok and ok_data
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
