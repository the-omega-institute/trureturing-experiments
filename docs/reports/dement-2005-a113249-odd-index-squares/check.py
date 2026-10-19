"""OEIS A113249 family conjecture: a(m, 2n+1) is a perfect square for all m, n.

a(m,0)=-1, a(m,1)=4, a(m,2)=-13+6(m-1)+3(m-1)^2, a(m,3)=(-8+m^2)^2,
a(m,n)=m^4 a(m,n-4)+(2m)^2 a(m,n-3)-4 a(m,n-1).
Checks (exact integers):
  1. the m=3 terms reproduce the A113249 data line;
  2. a(m,2r+1) = b(m,r)^2 with b(m,0)=2, b(m,1)=m^2-8, b(m,r+2)=-4 b(m,r+1)-m^2 b(m,r), for 0<=m<=60, r<60;
  3. the closed form a(m,n) = m^(n+1)(1-(-1)^n)/4 + c(m,n+1)/2 with c(m,0)=1, c(m,1)=-2,
     c(m,k+1)=-4 c(m,k)-m^2 c(m,k-1), and c(m,2k)=2 c(m,k)^2 - m^(2k), for the same range.
"""
from math import isqrt

def seq(m, N):
    a = [-1, 4, -13 + 6 * (m - 1) + 3 * (m - 1) ** 2, (-8 + m * m) ** 2]
    while len(a) < N:
        n = len(a)
        a.append(m ** 4 * a[n - 4] + (2 * m) ** 2 * a[n - 3] - 4 * a[n - 1])
    return a

DATA = [-1, 4, 11, 1, 59, 484, -1009, 6241, -2761, 13924, 87251, 57121, 49139, 4072324]
ok_data = seq(3, len(DATA)) == DATA
ok_sq = ok_closed = ok_double = True
for m in range(0, 61):
    a = seq(m, 122)
    b = [2, m * m - 8]
    while len(b) < 61:
        b.append(-4 * b[-1] - m * m * b[-2])
    c = [1, -2]
    while len(c) < 124:
        c.append(-4 * c[-1] - m * m * c[-2])
    for r in range(60):
        x = a[2 * r + 1]
        if x < 0 or isqrt(x) ** 2 != x or x != b[r] ** 2:
            ok_sq = False
    for n in range(122):
        num = m ** (n + 1) * (1 - (-1) ** n) + 2 * c[n + 1]
        if num % 4 != 0 or num // 4 != a[n]:
            ok_closed = False
    for k in range(61):
        if c[2 * k] != 2 * c[k] ** 2 - m ** (2 * k):
            ok_double = False
print("m=3 reproduces the A113249 data line:", ok_data)
print("a(m,2r+1)=b(m,r)^2 for 0<=m<=60, r<60:", ok_sq)
print("closed form via c for 0<=m<=60, n<122:", ok_closed)
print("c(m,2k)=2c(m,k)^2-m^(2k) for 0<=m<=60, k<=60:", ok_double)
ok = ok_data and ok_sq and ok_closed and ok_double
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
