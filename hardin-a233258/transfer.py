"""Column transfer of the A233258 three-row antidiagonal strip (requires sympy).

Labels are complementary-pair classes 0..2; consecutive columns c, d are compatible
when c[r] != d[r] for every row r and d[r] != c[r+1] (antidiagonal bond).
"""
import itertools
import sympy as sp

R = 3
cols = list(itertools.product(range(3), repeat=R))
def ok(c, d):
    return all(c[i] != d[i] for i in range(R)) and all(d[i] != c[i + 1] for i in range(R - 1))
T = sp.Matrix([[1 if ok(c, d) else 0 for d in cols] for c in cols])
lam = sp.symbols('lambda')
chi = sp.factor(T.charpoly(lam).as_expr())
assert sp.expand(chi - lam**15 * (lam - 2) * (lam - 1)**3 * (lam + 1)**2
                 * (lam**2 - 5*lam + 2) * (lam**2 + 4*lam + 2)**2) == 0

def proper(n):
    v = {c: 1 for c in cols}
    for _ in range(n - 1):
        v = {d: sum(v[c] for c in cols if ok(c, d)) for d in cols}
    return sum(v.values())
L = {n: proper(n) for n in range(1, 40)}
fails = [n for n in range(4, 40) if L[n] + 6 != 6*(L[n-1] + 6) - 7*(L[n-2] + 6) + 2*(L[n-3] + 6)]
assert fails == [4, 5]
x = sp.symbols('x')
assert sp.expand((x - 8) * (x**2 - 40*x + 128) - (x**3 - 48*x**2 + 448*x - 1024)) == 0
print('charpoly', chi)
print('L(n)+6 satisfies y^3-6y^2+7y-2 exactly for n >= 6; fails at', fails)
