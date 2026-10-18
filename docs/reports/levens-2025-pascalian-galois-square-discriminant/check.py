"""Levens, arXiv:2511.03082v1, Conjecture 6.4: "For n = 2m + 1, the Galois group of P_n(z) is the mth hyperoctahedral
group", where P_n(z) = sum_{k=0}^{n} C(n, floor((n-k)/2)) z^(n-k). For odd n = 2m+1, P_n(z) = (1+z) q(z^2) with
q(x) = sum_{j=0}^{m} C(2m+1, j) x^j (Theorem 6.2).

Checks (exact integer arithmetic, sympy):
- P_49(z) = (1+z) q(z^2) for m = 24;
- the differential identity (1+x) q'(x) - (2m+1) q(x) = -(m+1) c x^m, c = C(2m+1, m), for 1 <= m <= 30;
- Disc(q) = (-1)^(m(m-1)/2) (2m+1) (m+1)^(m-1) c^(m-2) for 2 <= m <= 30;
- Disc(q) for m = 24 is the square (7 * 5^23 * c^11)^2;
- among 2 <= m <= 40, Disc(q) is a square exactly for m = 24;
- the Pell family: m = (a^2-1)/2 with a^2 - 2 b^2 = -1 gives 2m+1 = a^2, m+1 = b^2 and m = 0 (mod 4), hence a square
  discriminant by the closed form (listed, not computed directly beyond m = 24).
Requires sympy.
"""
import sympy as sp
from math import comb, isqrt
x, z = sp.symbols('x z')
def q(m): return sum(comb(2*m+1, j)*x**j for j in range(m+1))
ok = True
m = 24; n = 2*m+1
P = sum(comb(n, (n-k)//2)*z**(n-k) for k in range(n+1))
ok &= sp.expand(P - (1+z)*q(m).subs(x, z**2)) == 0
print("P_49 = (1+z) q(z^2):", ok)
for mm in range(1, 31):
    c = comb(2*mm+1, mm)
    ok &= sp.expand((1+x)*sp.diff(q(mm), x) - (2*mm+1)*q(mm) + (mm+1)*c*x**mm) == 0
print("differential identity for m = 1..30:", ok)
form_ok = True
for mm in range(2, 31):
    c = comb(2*mm+1, mm)
    D = int(sp.discriminant(q(mm), x))
    form_ok &= D == (-1)**(mm*(mm-1)//2)*(2*mm+1)*(mm+1)**(mm-1)*c**(mm-2)
print("discriminant closed form for m = 2..30:", form_ok); ok &= form_ok
c = comb(49, 24); D24 = int(sp.discriminant(q(24), x)); r = 7*5**23*c**11
ok &= D24 == r*r
print("Disc(q_24) = (7*5^23*c^11)^2:", D24 == r*r)
squares = [mm for mm in range(2, 41) if (lambda D: D > 0 and isqrt(D)**2 == D)(int(sp.discriminant(q(mm), x)))]
print("square discriminants for m = 2..40:", squares); ok &= squares == [24]
a, b, fam = 1, 1, []
for _ in range(4):
    a, b = 3*a + 4*b, 2*a + 3*b; fam.append(((a*a-1)//2, a, b))
print("Pell family (m, sqrt(2m+1), sqrt(m+1)):", fam, " all m = 0 mod 4:", all(f[0] % 4 == 0 for f in fam))
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
