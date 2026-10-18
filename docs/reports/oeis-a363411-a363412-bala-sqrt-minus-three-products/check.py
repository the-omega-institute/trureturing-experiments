"""OEIS A363411, A363412 (P. Bala): write P_n = prod_{k=0..n} (1 + k*sqrt(-3)) = A_n + B_n*sqrt(-3).
A363411(n) = A_n (real part), A363412(n) = B_n (imaginary part divided by sqrt(3)).

Checked statements:
- A363412-1: 3 does not divide B_{3n+1} for all n.
- A363412-2: v_3(B_{3n}) = v_3(B_{3n-1}) for all n >= 1.
- A363411-1: v_5(A_{n+2}) = A079998(n) = [5 | n] for all n >= 0.
- A363411-2: the type 1 primes (primes dividing no A_n) up to 101 are exactly 3, 11, 17, 23, 29, 41, 47, 53, 59, 83, 89, 101.

Proof ingredients checked here:
- Lemma: for every m >= 1 and every odd j, T_m = m(m+1)/2 divides e_j(1, ..., m); hence
  B_m = sum_i e_{2i+1}(-3)^i = T_m (1 + 3K_m) and v_3(B_m) = v_3(T_m).
- For an inert prime p (here p = 2 mod 3), prod_{k=0..p-1} (1 + k t) = 2 (mod p) in Z[t]/(t^2+3),
  so P_{n+p} = 2 P_n (mod p) and A_n = 0 (mod p) has a solution iff it has one with n < p.
- prod_{k=0..24} (1 + k t) = 7 (mod 25), so P_{n+25} = 7 P_n (mod 25), and v_5(A_n) depends on n mod 25.
- 3 never divides A_n (A_n = 1 mod 3); a split prime p (p = 1 mod 6) and 2 divide some A_n (exhibited).

The script verifies the lemma for m <= L, the four statements for n <= N with exact integers, the
mod-p and mod-25 products, the finite residue checks, and the published data.
"""
import sys
from sympy import primerange
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
L = int(sys.argv[2]) if len(sys.argv) > 2 else 300
AD = [1,1,-5,-32,112,2212,-5348,-292880,276976,64180144,60400144]
BD = [0,1,3,-12,-140,420,13692,-23744,-2366784,126000,641927440,1306329024]
A, B = [1], [0]
for k in range(1, N + 1):
    a, b = A[-1], B[-1]
    A.append(a - 3 * k * b); B.append(b + k * a)
def v(x, p):
    if x == 0: return None
    c = 0
    while x % p == 0: x //= p; c += 1
    return c
ok_data = A[:len(AD)] == AD and B[:len(BD)] == BD
def esym(m):
    e = [1] + [0] * m
    for k in range(1, m + 1):
        for j in range(k, 0, -1): e[j] += k * e[j - 1]
    return e
ok_lemma = all(esym(m)[j] % (m * (m + 1) // 2) == 0 for m in range(1, L + 1) for j in range(1, m + 1, 2))
ok_b = all(B[m] % (3 * (m * (m + 1) // 2)) == m * (m + 1) // 2 % (3 * (m * (m + 1) // 2)) for m in range(1, N + 1))
c1 = all(B[3 * n + 1] % 3 != 0 for n in range(0, (N - 1) // 3 + 1) if 3 * n + 1 <= N)
c2 = all(v(B[3 * n], 3) == v(B[3 * n - 1], 3) for n in range(1, N // 3 + 1))
c3 = all(v(A[n + 2], 5) == (1 if n % 5 == 0 else 0) for n in range(0, N - 1))
def mul(x, y, mod):
    return ((x[0] * y[0] - 3 * x[1] * y[1]) % mod, (x[0] * y[1] + x[1] * y[0]) % mod)
def full(mod, length):
    X = (1, 0)
    for k in range(length): X = mul(X, (1, k), mod)
    return X
inert_ok = all(full(p, p) == (2, 0) for p in primerange(5, 102) if p % 3 == 2)
q25 = full(25, 25)
res25 = True
X = (1, 0); P = []
for k in range(25): X = mul(X, (1, k), 25); P.append(X)
for r in range(25):
    Y = P[r]
    for j in range(4):
        z5, z25 = Y[0] % 5 == 0, Y[0] % 25 == 0
        if z5 != (r % 5 == 2) or z25: res25 = False
        Y = mul(Y, q25, 25)
listed = [3, 11, 17, 23, 29, 41, 47, 53, 59, 83, 89, 101]
type1 = []
for p in primerange(2, 102):
    if p == 3:
        never = all(a % 3 == 1 for a in A)
    elif p % 3 == 2 and p != 2:
        never = all(a % p != 0 for a in A[:p])
    else:
        never = not any(a % p == 0 for a in A[:2 * p + 2])
    if never: type1.append(p)
c4 = type1 == listed and all((p in listed) == all(a % p != 0 for a in A) for p in primerange(2, 102))
print(f"data: {ok_data}; lemma T_m | e_odd for m<={L}: {ok_lemma}; B_m = T_m mod 3T_m for m<={N}: {ok_b}")
print(f"A363412-1 (3 does not divide B_(3n+1)): {c1}; A363412-2 (v_3 equality): {c2}  [n up to N/3]")
print(f"A363411-1 (v_5(A_(n+2)) = [5|n]) n<={N-2}: {c3}; product over a residue system mod 25 = {q25}; residue check: {res25}")
print(f"inert p<=101: product over F_p = 2: {inert_ok}; type 1 primes <= 101 = {type1}; A363411-2: {c4}")
good = ok_data and ok_lemma and ok_b and c1 and c2 and c3 and q25 == (7, 0) and res25 and inert_ok and c4
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
