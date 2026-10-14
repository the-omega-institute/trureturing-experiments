#!/usr/bin/env python3
"""Exact checks for Kraus-Simanek, arXiv:2309.16450v1 (J. Elasticity 156 (2024) 1085), Conjecture 2.1:
'For each N, there is an n so that among all n-gons with area 1, rho_N has no maximizer.'
rho_N(Omega) = min over polynomials p of degree <= N of the integral over Omega of |conj(z) - p(z)|^2 dA.

Construction: m = max(3, N+2), n = 2m, theta = pi/m, zeta = exp(2 pi i/m); R > 0 with R^2 > 1/(m sin theta cos theta),
r = 1/(m R sin theta); vertices v_{2j} = R zeta^j, v_{2j+1} = r e^{i theta} zeta^j. Fan triangles [0, v_j, v_{j+1}].
Exact (sympy, symbolic R) for N = 0..4:
1. area = 1;
2. the moments int z^k dA vanish for 1 <= k <= N+1 (checked to 1e-45 at 60 digits for R = 3/2, 2, 7), so rho_N = int |z|^2 dA (p = 0 optimal);
3. int |z|^2 dA = closed form that tends to infinity as R -> infinity (leading term (m sin(2 theta) / 12) * ... R^2-type growth).
Prints ALL_OK.
"""
import sympy as sp
R = sp.symbols('R', positive=True)
s, t = sp.symbols('s t', nonnegative=True)
ok = True
def tri_integral(f, a, b):
    # integral over triangle [0, a, b] of f(z, conj z) dA = |Im(conj(a) b)| * int_{s,t>=0, s+t<=1} f(s a + t b) ds dt
    jac = sp.Abs(sp.im(sp.conjugate(a) * b))
    z = s * a + t * b
    g = sp.expand(f(z, sp.conjugate(z)))
    inner = sp.integrate(g, (t, 0, 1 - s))
    return sp.simplify(jac * sp.integrate(inner, (s, 0, 1)))
for N in range(0, 5):
    m = max(3, N + 2); n = 2 * m
    th = sp.pi / m; zeta = sp.exp(2 * sp.pi * sp.I / m)
    r = 1 / (m * R * sp.sin(th))
    V = []
    for j in range(m):
        V.append(R * zeta ** j); V.append(r * sp.exp(sp.I * th) * zeta ** j)
    edges = [(V[k], V[(k + 1) % n]) for k in range(n)]
    area = sp.simplify(sum(tri_integral(lambda z, zb: 1, a, b) for a, b in edges))
    if sp.simplify(area - 1) != 0: print('AREA', N, area); ok = False
    for k in range(1, N + 2):
        Mk = sum(tri_integral(lambda z, zb, k=k: z ** k, a, b) for a, b in edges)
        # exact zero is decided numerically at 60 digits for three values of R (nested radicals of cos(pi/m) defeat simplify)
        vals = [abs(sp.N(Mk.subs(R, Rv), 60)) for Rv in (sp.Rational(3, 2), 2, 7)]
        if max(vals) > sp.Float('1e-45'): print('MOMENT', N, k, vals); ok = False
    I2 = sp.simplify(sp.expand_complex(sum(tri_integral(lambda z, zb: z * zb, a, b) for a, b in edges)))
    lim = sp.limit(I2, R, sp.oo)
    print('N=%d n=%d  int|z|^2 = %s   limit R->oo: %s' % (N, n, sp.simplify(I2), lim))
    if lim != sp.oo: ok = False
print('ALL_OK' if ok else 'FAIL')
