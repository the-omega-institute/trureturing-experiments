#!/usr/bin/env python3
"""Side entanglement helps to discriminate channels whose Kraus operators are all positive definite.

Source: M. Rexiti, L. Memarzadeh, S. Mancini, "Discrimination of dephasing channels", arXiv:2201.00388v2,
J. Phys. A 55, 245301 (2022), Sec. 7: "we are led to conjecture that side entanglement is useful in channel
discrimination if the Kraus operators of the channel admit a polar decomposition with non-trivial unitaries
(i.e., not identities). ... In contrast, if all Kraus operators have simply the identity as unitary part in the
polar decomposition, there is no benefit from side entanglement."

Family (a > sqrt 3): K_s(a) = (a I + s_x X + s_y Y + s_z Z) / (2 sqrt(a^2 + 3)), s in the tetrahedron
T = {(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)}. Exact checks (sympy) for a = 2 and a = 3:
  * every K_s is Hermitian positive definite (polar unitary factor = I), sum K_s^2 = I, the four are independent;
  * Phi_a(A) = lam_a A + (1 - lam_a) Tr(A) I / 2, lam_a = (a^2 - 1)/(a^2 + 3): lam_2 = 3/7, lam_3 = 2/3;
  * unassisted: (Phi_3 - Phi_2)(rho) = (5/21)(rho - I/2) for every density rho (checked symbolically), whose
    trace norm is at most 5/21, so every unassisted strategy succeeds with probability <= 1/2 + 5/84 = 47/84;
  * assisted with the Bell state P: (id x Phi_3 - id x Phi_2)(P) has eigenvalues 5/28, -5/84 (x3), and the
    measurement {P, I - P} succeeds with probability 1/2 + 5/56 = 33/56 > 47/84.
Prints ALL_OK.
"""
import sys
import sympy as sp

def fail(m):
    print("FAIL", m); sys.exit(1)

I2 = sp.eye(2); X = sp.Matrix([[0, 1], [1, 0]]); Y = sp.Matrix([[0, -sp.I], [sp.I, 0]]); Z = sp.Matrix([[1, 0], [0, -1]])
T = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
def kraus(a):
    return [(a * I2 + s[0] * X + s[1] * Y + s[2] * Z) / (2 * sp.sqrt(a**2 + 3)) for s in T]
def channel(a, A):
    return sp.simplify(sum((K * A * K.H for K in kraus(a)), sp.zeros(2, 2)))
lam = {}
for a in (2, 3):
    Ks = kraus(a)
    for K in Ks:
        if sp.simplify(K - K.H) != sp.zeros(2, 2): fail("not Hermitian")
        ev = [sp.nsimplify(sp.simplify(e)) for e in K.eigenvals()]
        if not all(sp.simplify(e) > 0 for e in ev): fail(f"not positive definite a={a}")
    if sp.simplify(sum((K * K for K in Ks), sp.zeros(2, 2)) - I2) != sp.zeros(2, 2): fail("not trace preserving")
    coeff = sp.Matrix([[a, *s] for s in T])
    if coeff.det() == 0: fail("dependent Kraus operators")
    p, q, r, t = sp.symbols('p q r t')
    A = sp.Matrix([[p, q], [r, t]])
    la = sp.Rational(a**2 - 1, a**2 + 3); lam[a] = la
    if sp.simplify(channel(a, A) - (la * A + (1 - la) * A.trace() * I2 / 2)) != sp.zeros(2, 2): fail("channel action")
if (lam[2], lam[3]) != (sp.Rational(3, 7), sp.Rational(2, 3)): fail("lambda values")
delta = lam[3] - lam[2]
if delta != sp.Rational(5, 21): fail("delta")
# unassisted bound: difference is delta (rho - I/2); its eigenvalues are +/- delta * |2 r - 1| / 2 <= delta / 2
pun = sp.Rational(1, 2) + delta / 4
if pun != sp.Rational(47, 84): fail("unassisted optimum")
# assisted: Bell projector
e00 = sp.Matrix([1, 0, 0, 1]) / sp.sqrt(2); P = e00 * e00.H
def id_kron(a):
    out = sp.zeros(4, 4)
    for K in kraus(a):
        KK = sp.kronecker_product(I2, K)
        out += KK * P * KK.H
    return sp.simplify(out)
Dif = sp.simplify(id_kron(3) - id_kron(2))
ev = Dif.eigenvals()
if ev != {sp.Rational(5, 28): 1, sp.Rational(-5, 84): 3}: fail(f"assisted eigenvalues {ev}")
pbell = sp.Rational(1, 2) + sp.simplify((P * Dif).trace()) / 2
if pbell != sp.Rational(33, 56) or not pbell > pun: fail(f"assisted success {pbell}")
print(f"lam_2 = {lam[2]}, lam_3 = {lam[3]}, delta = {delta}")
print(f"unassisted optimum = {pun}; Bell-assisted success = {pbell} = 1/2 + 5/56 > 47/84")
print("ALL_OK")
