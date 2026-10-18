"""OEIS A377106 (Hanna): G.f. A(x)^3 = A(x^3 + 9 x A(x)^3), A(0) = 0, A'(0) = 1. Conjecture: a(n) is odd iff n = 2^k.

Route: the reversion R = A^(-1) satisfies R(y^3) = R(y)^3 + 9 y^3 R(y). Its coefficient recursion is
3 r_n = [y^(n+2)] (R(y^3) - R_<n(y)^3 - 9 y^3 R_<n(y)), and the bracket is divisible by 3 because
R_<n(y)^3 = R_<n(y^3) (mod 3) for integer coefficients. So R is integral, R = y + y^2 (mod 2), and A is the
compositional inverse, A = sum_k x^(2^k) (mod 2).
Checks to order N: every division by 3 is exact; R mod 2 = y + y^2; A (inverse of R) agrees with the OEIS
data; A satisfies the defining equation; a(n) is odd exactly at powers of 2.
"""
N = 150

def mul(a, b, n):
    c = [0] * (n + 1)
    for i, x in enumerate(a[:n + 1]):
        if x:
            for j, y in enumerate(b[:n + 1 - i]):
                if y: c[i + j] += x * y
    return c

def compose(f, g, n):
    # f(g(x)) for g(0) = 0, coefficients up to x^n
    res = [0] * (n + 1); P = [1] + [0] * n
    for k in range(1, n + 1):
        P = mul(P, g, n)
        if f[k]:
            for i in range(n + 1): res[i] += f[k] * P[i]
    return res

ok = True
r = [0] * (N + 3); r[1] = 1; exact = True
for n in range(2, N + 1):
    Rt = r[:n] + [0] * (N + 3 - n)
    R3 = mul(mul(Rt, Rt, n + 2), Rt, n + 2)
    lhs = r[(n + 2) // 3] if (n + 2) % 3 == 0 else 0
    rest = lhs - R3[n + 2] - 9 * r[n - 1]
    exact &= rest % 3 == 0
    r[n] = rest // 3
print('R integral to order', N, '(every division by 3 exact):', exact); ok &= exact
mod2 = [x % 2 for x in r[1:N + 1]]
print('R mod 2 = y + y^2:', mod2 == [1, 1] + [0] * (N - 2)); ok &= mod2 == [1, 1] + [0] * (N - 2)
A = [0] * (N + 1); A[1] = 1
for n in range(2, N + 1):
    A[n] = -compose(r, A, n)[n]
oeis = [1, 3, 18, 127, 966, 7686, 63068, 529503, 4526262, 39262658, 344789172, 3059733222]
print('first 12 terms match the OEIS data:', A[1:13] == oeis); ok &= A[1:13] == oeis
A3 = mul(mul(A, A, N), A, N)
inner = [0] * (N + 1); inner[3] += 1
for i, v in enumerate(mul([0, 9] + [0] * (N - 1), A3, N)): inner[i] += v
rhs = compose(A, inner, N)
print('A^3 = A(x^3 + 9 x A^3) to order', N, ':', A3 == rhs); ok &= A3 == rhs
odd = [n for n in range(1, N + 1) if A[n] % 2]
print('odd a(n), n <=', N, ':', odd); ok &= odd == [1, 2, 4, 8, 16, 32, 64, 128]
print('RESULT', 'PASS' if ok else 'FAIL')
