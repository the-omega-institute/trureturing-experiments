#!/usr/bin/env python3
"""A purely non-unital, eternally non-Markovian qubit dynamics.

Source: S. Utagi, S. Banerjee, R. Srikanth, "On the eternal non-Markovianity of non-unital quantum channels",
arXiv:2203.10625v2, Int. J. Quantum Inf. 22, 2350039 (2024): Definition 1 (non-Markovianity originating solely from
the translation tau(t) of r -> M(t) r + tau(t)) and Sec. IV-V: "we conjecture that for a non-unital channel with no
such unital component (in the above sense) to be ENM, is impossible."

Dynamics: D_t = 1 + t^2, P_t = [[1, t], [t, t^2]] / D_t, Phi_t(X) = (X + t^2 Tr(X) P_t) / D_t.
Exact checks (sympy):
  1. Kraus operators I/D, tI/D, (t/D)[[1,0],[t,0]], (t/D)[[0,1],[0,t]] are complete and give Phi_t; Phi_0 = id;
     Phi_t(I) != I for t > 0.
  2. The unital part Delta_t(X) = Phi_t(X) - Tr(X)(Phi_t(I) - I)/2 = (X + t^2 Tr(X) I/2)/D_t is CP-divisible:
     Delta_t = W_{t,s} o Delta_s with W_{t,s} = eta id + (1 - eta) Tr(.) I/2, eta = D_s/D_t in (0, 1].
  3. The unique intermediate map V_{t,s} = Phi_t o Phi_s^{-1} is not positive for 0 < s < t:
     (-t, 1) V_{t,s}(P_t) (-t, 1)^T = -s^2 (t - s)^2 / (D_s D_t) < 0.
  4. The generator L_t = dPhi_t/dt o Phi_t^{-1} has a Choi matrix with a negative eigenvalue on the complement of
     the maximally entangled vector (the minimum canonical rate), equal to 40 digits to
     (t/D_t)(1 - sqrt(1 + t^2/D_t^2)) at t in {1/10, 1/2, 1, 2, 5}.
Prints ALL_OK.
"""
import sys
import sympy as sp
import mpmath as mp

def fail(m):
    print("FAIL", m); sys.exit(1)

t, s = sp.symbols('t s', positive=True)
I2 = sp.eye(2)
def E(i, j):
    m = sp.zeros(2, 2); m[i, j] = 1; return m
def P(x): return sp.Matrix([[1, x], [x, x**2]]) / (1 + x**2)
def Phi(x, X): return (X + x**2 * X.trace() * P(x)) / (1 + x**2)
a, b, c, d = sp.symbols('a b c d'); X = sp.Matrix([[a, b], [c, d]])
D = 1 + t**2
K = [I2 / D, t * I2 / D, (t / D) * sp.Matrix([[1, 0], [t, 0]]), (t / D) * sp.Matrix([[0, 1], [0, t]])]
if sp.simplify(sum((k.H * k for k in K), sp.zeros(2, 2)) - I2) != sp.zeros(2, 2): fail("Kraus completeness")
if sp.simplify(sum((k * X * k.H for k in K), sp.zeros(2, 2)) - Phi(t, X)) != sp.zeros(2, 2): fail("Kraus action")
if sp.simplify(Phi(0, X) - X) != sp.zeros(2, 2): fail("Phi_0")
if sp.simplify(Phi(t, I2) - I2) == sp.zeros(2, 2): fail("unital")
def Delta(x, Y): return Phi(x, Y) - Y.trace() * (Phi(x, I2) - I2) / 2
if sp.simplify(Delta(t, X) - (X + t**2 * X.trace() * I2 / 2) / D) != sp.zeros(2, 2): fail("unital part")
eta = (1 + s**2) / (1 + t**2)
W = lambda Y: eta * Y + (1 - eta) * Y.trace() * I2 / 2
if sp.simplify(W(Delta(s, X)) - Delta(t, X)) != sp.zeros(2, 2): fail("unital CP-divisibility")
# V_{t,s} = Phi_t o Phi_s^{-1}, Phi_s^{-1}(Y) = D_s Y - s^2 Tr(Y) P_s
Phis_inv = lambda Y: (1 + s**2) * Y - s**2 * Y.trace() * P(s)
if sp.simplify(Phi(s, Phis_inv(X)) - X) != sp.zeros(2, 2): fail("inverse")
Vout = Phi(t, Phis_inv(P(t)))
v = sp.Matrix([-t, 1])
q = sp.simplify((v.T * Vout * v)[0])
if sp.simplify(q + s**2 * (t - s)**2 / ((1 + s**2) * (1 + t**2))) != 0: fail("negative quadratic form")
print("1-3: Kraus/CPTP, unital part CP-divisible, V_{t,s}(P_t) has form -s^2(t-s)^2/(D_s D_t) (exact)")
# 4: canonical rates
mp.mp.dps = 40
basis = [(0, 0), (0, 1), (1, 0), (1, 1)]
def sup(f):
    return sp.Matrix.hstack(*[sp.Matrix(list(f(E(i, j)))) for (i, j) in basis])
St = sup(lambda Y: Phi(t, Y))
L = sp.simplify(sp.diff(St, t) * St.inv())
J = sp.zeros(4, 4)
for k, (i, j) in enumerate(basis):
    col = L[:, k]; J += sp.kronecker_product(E(i, j), sp.Matrix([[col[0], col[1]], [col[2], col[3]]]))
Om = sp.Matrix([1, 0, 0, 1]) / sp.sqrt(2); Pp = sp.eye(4) - Om * Om.T
for tv in [sp.Rational(1, 10), sp.Rational(1, 2), 1, 2, 5]:
    Jn = mp.matrix((Pp * J * Pp).subs(t, tv).evalf(50).tolist())
    Jn = (Jn + Jn.H) / 2
    ev = mp.eigsy(mp.matrix([[mp.re(Jn[i, j]) for j in range(4)] for i in range(4)]))[0] if all(abs(mp.im(Jn[i, j])) < mp.mpf(10)**-30 for i in range(4) for j in range(4)) else None
    if ev is None:
        ev = mp.eighe(Jn)[0]
    emin = min(ev)
    Dv = 1 + tv**2
    gmin = mp.mpf(tv) / Dv * (1 - mp.sqrt(1 + mp.mpf(tv)**2 / Dv**2))
    if not (emin < 0 and abs(emin - gmin) < mp.mpf(10)**-25): fail(f"rate t={tv}: {emin} vs {gmin}")
    print(f"4: t={tv}: minimum canonical rate {mp.nstr(emin, 15)} = (t/D)(1 - sqrt(1 + t^2/D^2)) < 0")
print("ALL_OK")
