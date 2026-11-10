"""Column transfer of the A233259 four-row antidiagonal strip (requires sympy).

Label columns Fin 4 -> Fin 3; consecutive columns c, d are compatible when
c[r] != d[r] for every row and d[r] != c[r+1] (antidiagonal bond).
"""
import itertools
import sympy as sp

R = 4
cols = list(itertools.product(range(3), repeat=R))
def ok(c, d):
    return all(c[i] != d[i] for i in range(R)) and all(d[i] != c[i + 1] for i in range(R - 1))
lam, y, x = sp.symbols('lambda y x')
T = sp.Matrix([[1 if ok(c, d) else 0 for d in cols] for c in cols])
chi = T.charpoly(lam).as_expr()
expected = (lam**65 * (lam - 2) * (lam + 1)**2 * (lam**2 - 6)**2 * (lam**2 - 5*lam + 5)
            * (lam**2 + 7*lam + 8)**2 * (lam**3 - 9*lam**2 + 15*lam - 6))
assert sp.expand(chi - expected) == 0

def proper(n):
    v = {c: 1 for c in cols}
    for _ in range(n - 1):
        v = {d: sum(v[c] for c in cols if ok(c, d)) for d in cols}
    return sum(v.values())
L = {n: proper(n) for n in range(1, 30)}
rel = lambda n: L[n] + 6 - (10*(L[n-1] + 6) - 24*(L[n-2] + 6) + 21*(L[n-3] + 6) - 6*(L[n-4] + 6))
fails = [n for n in range(5, 30) if rel(n) != 0]
assert fails == [5, 6, 7, 8]
assert sp.expand((y - 1) * (y**3 - 9*y**2 + 15*y - 6) - (y**4 - 10*y**3 + 24*y**2 - 21*y + 6)) == 0
assert sp.expand((x**4 - 160*x**3 + 6144*x**2 - 86016*x + 393216) - 16**4 * (y**4 - 10*y**3 + 24*y**2 - 21*y + 6).subs(y, x / 16)) == 0
print('charpoly', sp.factor(chi))
print('L(n)+6 relation (y-1)(y^3-9y^2+15y-6) holds exactly for n >= 9; fails at', fails)
