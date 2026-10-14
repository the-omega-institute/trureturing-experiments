#!/usr/bin/env python3
"""Global optimality of the aligned general-dyne measurement for a monitored thermal OPO.

Source: K. H. Kua, A. Serafini, M. G. Genoni, "Daemonic ergotropy of Gaussian quantum states and the role of
measurement-induced purification via general-dyne detection", arXiv:2506.22288v2, Sec. VI, after Eq. (66):
"We have also numerically verified that non-zero phases theta_m lead to worse result, and we can thus conjecture
that this fully characterizes the optimal steady-state general-dyne monitoring."

Dimensionless steady state of Eq. (55) for the OPO (57) with A = diag(-k/2-chi, -k/2+chi), D = k nu I and a
passive coupling (r = 2 chi / k in [0,1), nu = 2 n_th + 1 >= 1, M = sigma_m the measurement covariance):
    0 = nu I - (D0 S + S D0)/2 - (S - nu I) Q (S - nu I),  D0 = diag(1+r, 1-r),  Q = (nu I + M)^(-1).
Claimed optimum: M* = diag((1-r)/(1+r), (1+r)/(1-r)) (theta = 0, z = (1-r)/(1+r), pure), with
    sqrt(det S) >= v* = [(nu-1) c + sqrt((nu-1)^2 c^2 + 4 nu)] / 2,   c = sqrt(1 - r^2).
Part 1 (sympy, exact): the two determinant identities of the proof sketch.
Part 2 (float64): integrate the Riccati flow to its stabilizing steady state for 3000 random physical M
  (rotated, squeezed, noisy, near-homodyne) over a grid of (r, nu); check sqrt(det S) >= v* - 1e-9 and equality at M*.
Prints ALL_OK.
"""
import sys
import numpy as np
import sympy as sp

def fail(m):
    print("FAIL", m); sys.exit(1)

# Part 1
x, y, h, nu, r = sp.symbols('x y h nu r', real=True)
a, b = 1 + r, 1 - r
S = sp.Matrix([[x, h], [h, y]]); I2 = sp.eye(2); D0 = sp.diag(a, b)
V = S - nu * I2; R = nu * I2 - (D0 * S + S * D0) / 2
d = x * y - h**2; T = a * x + b * y
lhs1 = sp.expand((V * V - nu * R).det() - R.det())
rhs1 = sp.expand((d - 1) * (d + nu**2 - nu * T) - (nu**2 - 1) * r**2 * x * y)
if sp.simplify(lhs1 - rhs1) != 0: fail("identity 1")
v, c, eps = sp.symbols('v c eps', positive=True)
Fp = v**2 - (nu + 1) * c * v + nu; Fm = v**2 - (nu - 1) * c * v - nu
# substitute d = v^2, c^2 = a b, and identity 1 with eps = (det M - 1) det R
lhs2 = sp.expand(Fp * Fm)
rhs2 = sp.expand((nu**2 - 1) * r**2 * h**2 + nu * (v**2 - 1) * (T - 2 * c * v) + eps)
# eps = (v^2-1)(v^2+nu^2-nu T) - (nu^2-1) r^2 x y  from identity 1 at d = v^2, and x y = v^2 + h^2
eps_expr = (v**2 - 1) * (v**2 + nu**2 - nu * T) - (nu**2 - 1) * r**2 * (v**2 + h**2)
diff = sp.expand((lhs2 - rhs2.subs(eps, eps_expr)).subs(c**2, a * b))
diff = sp.expand(diff.subs(c**3, a * b * c).subs(c**2, a * b))
if sp.simplify(diff) != 0: fail("identity 2")
print("part1: det(V^2 - nu R) - det R = (d-1)(d+nu^2-nu T) - (nu^2-1) r^2 x y and "
      "F+(v) F-(v) = (nu^2-1) r^2 h^2 + nu (d-1)(T - 2 c v) + eps (exact)")

# Part 2
rng = np.random.default_rng(250622288)
def steady(rv, nuv, M):
    D0n = np.diag([1 + rv, 1 - rv]); Q = np.linalg.inv(nuv * np.eye(2) + M); Sn = nuv * np.eye(2)
    dt = 0.02
    for _ in range(200000):
        Vn = Sn - nuv * np.eye(2)
        F = nuv * np.eye(2) - (D0n @ Sn + Sn @ D0n) / 2 - Vn @ Q @ Vn
        Sn = Sn + dt * F; Sn = (Sn + Sn.T) / 2
        if np.abs(F).max() < 1e-13: break
    return Sn, np.abs(F).max()
count = 0
for rv in (0.1, 0.5, 0.9, 0.99):
    for nuv in (1.0, 1.5, 3.0, 10.0):
        cv = np.sqrt(1 - rv**2)
        vstar = ((nuv - 1) * cv + np.sqrt((nuv - 1)**2 * cv**2 + 4 * nuv)) / 2
        Mstar = np.diag([(1 - rv) / (1 + rv), (1 + rv) / (1 - rv)])
        Ss, res = steady(rv, nuv, Mstar)
        if res > 1e-9 or abs(np.sqrt(np.linalg.det(Ss)) - vstar) > 1e-7: fail(f"optimum r={rv} nu={nuv}")
        for _ in range(190):
            th = rng.uniform(0, np.pi); z = np.exp(rng.uniform(-6, 6)); nm = 1 + rng.exponential(1.0) * (rng.random() < 0.5)
            Rt = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
            M = nm * Rt @ np.diag([z, 1 / z]) @ Rt.T
            Sm, res = steady(rv, nuv, M)
            if res > 1e-8: fail(f"no convergence r={rv} nu={nuv}")
            if np.linalg.eigvalsh(Sm).min() <= 0: fail("non-PD steady state")
            if np.sqrt(np.linalg.det(Sm)) < vstar - 1e-9: fail(f"bound violated r={rv} nu={nuv} M={M}")
            count += 1
        print(f"part2: r={rv} nu={nuv}: v*={vstar:.10f} attained at M*; 190 random M respect the bound")
print(f"part2: {count} random measurements checked")
print("ALL_OK")
