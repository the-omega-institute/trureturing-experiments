"""Factor Hardin's order-7 polynomial for A233171 over the rationals (requires sympy)."""
import sympy as sp

x, y = sp.symbols('x y')
p = (x**7 - 480*x**6 + 60416*x**5 - 1966080*x**4 - 22020096*x**3
     + 1845493760*x**2 - 22548578304*x + 68719476736)
sextic = y**6 - 14*y**5 + 45*y**4 - 15*y**3 - 36*y**2 + 19*y - 2
assert sp.expand(p - 32**7 * ((y - 1) * sextic).subs(y, x / 32)) == 0
assert sp.Poly(sextic, y).is_irreducible
print('p(x) = 32^7 (y-1) s(y) with y = x/32, s irreducible over Q:', sextic)
