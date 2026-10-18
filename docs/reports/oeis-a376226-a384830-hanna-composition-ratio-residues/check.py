"""OEIS A376226, A384830 (Hanna): residues of the solutions of U(x^3 + a x U^3) = U * U(x^2 + b x U^2), U = x + ...
A376226: (a, b) = (6, 4); conjectures (C1) a(n) = 1 mod 3 iff n = 3^k, (C2) a(n) = 2 mod 3 iff 2n = 3^i + 3^j, i < j.
A384830: (a, b) = (-3, -2); (C.1) a(n) odd iff n = 2^k, (C.2), (C.3) as (C1), (C2).
Route (GPT Pro r386q, checked): the coefficient recursion has multiplier -1 (no division, unique over Z, F_2, F_3);
mod 2 the A384830 equation is solved by S = sum x^(2^k); mod 3 both equations become B(x^3) = B * B(x^2 + x B^2),
solved by B = inverse of y + y^2 + y^3, and B(t^2) = (sum t^(3^k))^2 gives the residues.
"""
def mul(a, b, n):
    c = [0] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j, y in enumerate(b[:n + 1 - i]):
                if y: c[i + j] += x * y
    return c
def power(a, k, n):
    r = [1] + [0] * n
    for _ in range(k): r = mul(r, a, n)
    return r
def compose(f, g, n):
    res = [0] * (n + 1); P = [1] + [0] * n
    for k in range(1, n + 1):
        P = mul(P, g, n)
        if k < len(f) and f[k]:
            for i in range(n + 1): res[i] += f[k] * P[i]
    return res
def E(U, al, be, M):
    in1 = [0] * (M + 1); in1[3] += 1
    for i, v in enumerate(mul([0, al] + [0] * (M - 1), power(U, 3, M), M)): in1[i] += v
    in2 = [0] * (M + 1); in2[2] += 1
    for i, v in enumerate(mul([0, be] + [0] * (M - 1), power(U, 2, M), M)): in2[i] += v
    L = compose(U + [0] * (M + 1), in1, M); R = mul(U, compose(U + [0] * (M + 1), in2, M), M)
    return [l - r for l, r in zip(L, R)]
def solve(al, be, N):
    U = [0] * (N + 1); U[1] = 1; mult_ok = True
    for m in range(2, N + 1):
        M = m + 2; Uc = U[:m] + [0] * (N + 1 - m)
        e0 = E(Uc, al, be, M)[M]; Uc[m] = 1; e1 = E(Uc, al, be, M)[M]
        mult_ok &= (e1 - e0) == -1
        U[m] = e0
    return U, mult_ok
def oeis(path, k):
    v = []
    for line in open(path):
        if line[:2] in ('%S', '%T', '%U'): v += [int(t) for t in line.split(None, 2)[2].strip().strip(',').split(',') if t]
    return v[:k]
N = 70
S2 = set((3 ** i + 3 ** j) // 2 for i in range(8) for j in range(i + 1, 8))
pow3 = lambda n: any(n == 3 ** k for k in range(20)); pow2 = lambda n: n & (n - 1) == 0
ok = True
F, m1 = solve(6, 4, N); H, m2 = solve(-3, -2, N)
print('recursion multiplier is -1 at every step:', m1 and m2); ok &= m1 and m2
for path, U in (('oeis-a376226.txt', F), ('oeis-a384830.txt', H)):
    ref = oeis(path, 12); m = U[1:13] == ref; print(f'{path}: first 12 terms match: {m}'); ok &= m
c1 = all((F[n] % 3 == 1) == pow3(n) and (F[n] % 3 == 2) == (n in S2) for n in range(1, N + 1))
c2 = all((H[n] % 2 == 1) == pow2(n) and (H[n] % 3 == 1) == pow3(n) and (H[n] % 3 == 2) == (n in S2) for n in range(1, N + 1))
print(f'A376226 (C1), (C2) for n <= {N}: {c1}'); print(f'A384830 (C.1)-(C.3) for n <= {N}: {c2}'); ok &= c1 and c2
# the explicit reductions: B = inverse of y + y^2 + y^3 over F_3 and S = sum x^(2^k) over F_2
h = [0, 1, 1, 1] + [0] * (N - 3); B = [0] * (N + 1); B[1] = 1
for m in range(2, N + 1):
    B[m] = (-compose(h, B, m)[m]) % 3
c3 = all(F[n] % 3 == B[n] % 3 == H[n] % 3 for n in range(1, N + 1))
c4 = all(H[n] % 2 == (1 if pow2(n) else 0) for n in range(1, N + 1))
print('F = H = inverse of y + y^2 + y^3 (mod 3):', c3, '  H = sum x^(2^k) (mod 2):', c4); ok &= c3 and c4
print('RESULT', 'PASS' if ok else 'FAIL')
