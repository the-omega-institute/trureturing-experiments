"""Factor Hardin's order-10 polynomial for A233172 over the rationals (requires sympy)."""
import sympy as sp

x, y = sp.symbols('x y')
c = [1216, -303104, -6553600, 6056574976, -287762808832, -18622978195456, 1482141674242048,
     -15199648742375424, -522417556774977536, 5764607523034234880]
p = x**10 - sum(ci * x**(9 - i) for i, ci in enumerate(c))
q = sp.expand(p.subs(x, 64 * y) / 64**10)
assert all(v.is_integer for v in sp.Poly(q, y).coeffs())
n = y**9 - 18*y**8 + 56*y**7 + 81*y**6 - 280*y**5 - 12*y**4 + 259*y**3 - 78*y**2 - 24*y + 5
assert sp.expand(q - (y - 1) * n) == 0 and sp.Poly(n, y).is_irreducible
print('p(x) = 64^10 (y-1) n(y) with y = x/64, n irreducible over Q:', n)
print('rational roots of p:', [r for r in sp.roots(sp.Poly(p, x), filter='Q')])
