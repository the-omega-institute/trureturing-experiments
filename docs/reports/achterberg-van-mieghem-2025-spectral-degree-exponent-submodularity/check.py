"""Check that the spectral degree exponent (arXiv:2502.01815v1, Definition 2.1) is not submodular on edge sets.

Vertex set {0,...,4}; E0 = {01, 02, 03, 14, 23}, e = 24, f = 34.  Part (a) gives exact characteristic polynomials and
rational certificates that compare lambda_1 with the third power mean M_3 of the degrees, so q(H0) > 3,
q(He) = q(Hf) < 3 and q(Hef) > 3.  Part (b) computes the exponents numerically, for the graphs and for their uniform
t-blowups (t = 2, 3).  Part (c) checks the induced-vertex-set reading and part (d) a supermodularity violation.
"""
from fractions import Fraction as Fr
import itertools
import numpy as np
import sympy as sp

z = sp.symbols('z')


def adj(n, edges):
    A = sp.zeros(n, n)
    for i, j in edges:
        A[i, j] = A[j, i] = 1
    return A


E0 = [(0, 1), (0, 2), (0, 3), (1, 4), (2, 3)]
e, f = (2, 4), (3, 4)
graphs = {'H0': E0, 'He': E0 + [e], 'Hf': E0 + [f], 'Hef': E0 + [e, f]}
ok = True
print("(a) exact certificates")
cert = {}
for name, E in graphs.items():
    A = adj(5, E)
    chi = sp.factor(A.charpoly(z).as_expr())
    deg = [sum(A.row(i)) for i in range(5)]
    m3cubed = sp.Rational(sum(d ** 3 for d in deg), 5)
    print(f"  {name}: degrees {sorted(deg, reverse=True)}, charpoly {chi}, M_3^3 = {m3cubed}")
    cert[name] = (A, sp.Poly(A.charpoly(z).as_expr(), z), deg, m3cubed)
# H0: chi0(11/5) < 0 and chi0 monic => root > 11/5; (11/5)^3 > M3^3 => lambda1 > M3 => q > 3
c0 = cert['H0'][1]
v = c0.eval(sp.Rational(11, 5)); w = sp.Rational(11, 5) ** 3 - cert['H0'][3]
print(f"  H0: chi(11/5) = {v} < 0; (11/5)^3 - M3^3 = {w} > 0 => q(H0) > 3"); ok &= v < 0 and w > 0
# He: largest root of a(z) = z^3 - 2z^2 - 2z + 2 lies below 249/100: a increasing on [2, oo), a(249/100) > 0
a = sp.Poly(z ** 3 - 2 * z ** 2 - 2 * z + 2, z)
assert sp.expand(cert['He'][1].as_expr() - z * (z + 2) * a.as_expr()) == 0
va = a.eval(sp.Rational(249, 100)); w = cert['He'][3] - sp.Rational(249, 100) ** 3
deriv_pos = sp.Poly(sp.diff(a.as_expr(), z), z).eval(2) > 0 and sp.diff(a.as_expr(), z, 2).subs(z, 2) > 0
print(f"  He: a(249/100) = {va} > 0, a' increasing and positive at 2; M3^3 - (249/100)^3 = {w} > 0 => q(He) < 3")
ok &= va > 0 and w > 0 and deriv_pos
# Hef: u(z) = z^3 - z^2 - 6z + 2, u(571/200) < 0 => root > 571/200; (571/200)^3 > M3^3 => q(Hef) > 3
u = sp.Poly(z ** 3 - z ** 2 - 6 * z + 2, z)
assert sp.expand(cert['Hef'][1].as_expr() - z * (z + 1) * u.as_expr()) == 0
vu = u.eval(sp.Rational(571, 200)); w = sp.Rational(571, 200) ** 3 - cert['Hef'][3]
print(f"  Hef: u(571/200) = {vu} < 0; (571/200)^3 - M3^3 = {w} > 0 => q(Hef) > 3"); ok &= vu < 0 and w > 0
print("  => q(He) + q(Hf) < 6 < q(H0) + q(Hef): submodularity fails")


def sde(Am):
    lam = max(np.linalg.eigvalsh(Am)); d = Am.sum(1)
    M = lambda q: (np.mean(d ** q)) ** (1 / q)
    lo, hi = 2.0, 2.0
    while M(hi) < lam:
        hi *= 2
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if M(mid) < lam else (lo, mid)
    return lam, (lo + hi) / 2


print("(b) numerical exponents")
for t in (1, 2, 3):
    vals = []
    for name, E in graphs.items():
        Eb = [(x * t + i, y * t + j) for x, y in E for i in range(t) for j in range(t)]
        Am = np.zeros((5 * t, 5 * t))
        for x, y in Eb:
            Am[x, y] = Am[y, x] = 1
        vals.append(sde(Am)[1])
    print(f"  t={t}: q(H0), q(He), q(Hf), q(Hef) = {', '.join('%.6f' % x for x in vals)}; "
          f"q(He)+q(Hf) - q(H0)-q(Hef) = {vals[1] + vals[2] - vals[0] - vals[3]:.6f}")
    ok &= vals[1] + vals[2] < vals[0] + vals[3]
print("(c) induced-vertex-set reading")
V = 'abcxyzw'; EG = ['ax', 'ay', 'az', 'aw', 'bx', 'by', 'bz', 'cx', 'cy', 'cw']


def induced(S):
    S = sorted(S, key=V.index); m = {c: i for i, c in enumerate(S)}
    Am = np.zeros((len(S), len(S)))
    for uv in EG:
        if uv[0] in m and uv[1] in m:
            Am[m[uv[0]], m[uv[1]]] = Am[m[uv[1]], m[uv[0]]] = 1
    return Am


S, T = set('abxyz'), set('acxyw')
qs = {k: sde(induced(X))[1] for k, X in (('S', S), ('T', T), ('S&T', S & T), ('S|T', S | T))}
print("  " + ", ".join(f"q(G[{k}]) = {v:.6f}" for k, v in qs.items()))
ok &= qs['S'] + qs['T'] < qs['S&T'] + qs['S|T']
print("(d) a supermodularity violation on 6 vertices")
E = [(0, 1), (4, 5), (3, 5), (2, 4), (1, 5), (3, 4), (1, 3)]; e2, f2 = (0, 5), (2, 5)
vals = []
for EE in (E, E + [e2], E + [f2], E + [e2, f2]):
    Am = np.zeros((6, 6))
    for x, y in EE:
        Am[x, y] = Am[y, x] = 1
    vals.append(sde(Am)[1])
print(f"  q(E), q(E+e), q(E+f), q(E+e+f) = {', '.join('%.6f' % x for x in vals)}; q(E+e)+q(E+f) > q(E)+q(E+e+f): {vals[1] + vals[2] > vals[0] + vals[3]}")
ok &= vals[1] + vals[2] > vals[0] + vals[3]
print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
