"""Factor the column recurrence polynomials of A233152, A233153, A233154 over Q and check Barker's g.f. numerators (requires sympy)."""
import sympy as sp

x = sp.symbols('x')
polys = {'A233152': x**2 - 15*x + 18, 'A233153': x**3 - 25*x**2 + 90*x - 81,
         'A233154': x**4 - 42*x**3 + 351*x**2 - 972*x + 810}
for k, p in polys.items():
    print(k, sp.factor_list(p), 'real roots:', [sp.N(r, 8) for r in sp.real_roots(p)])
assert sp.expand(polys['A233154'] - (x - 3) * (x**3 - 39*x**2 + 234*x - 270)) == 0
assert sp.Poly(x**3 - 39*x**2 + 234*x - 270, x).is_irreducible
assert sp.Poly(polys['A233153'], x).is_irreducible and sp.Poly(polys['A233152'], x).is_irreducible
