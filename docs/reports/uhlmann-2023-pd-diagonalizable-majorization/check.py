"""Uhlmann arXiv:2309.01302v1, Section 4: does diag(P D P^-1) majorize the spectrum diag(D) for every real
symmetric positive-definite P in every dimension? Exact (sympy) check of the 7x7 counterexample and family.
"""
import sys
import sympy as sp
P = sp.zeros(7, 7); P[0, 0] = P[1, 1] = 10; P[0, 1] = P[1, 0] = 1
for i in range(5):
    P[0, 2 + i] = P[2 + i, 0] = P[1, 2 + i] = P[2 + i, 1] = 2; P[2 + i, 2 + i] = 10
fails = []
if not (P == P.T and all(ev > 0 for ev in P.eigenvals())): fails.append("P not PD")
x, y = sp.symbols('x y'); z = sp.symbols('z0:5'); v = sp.Matrix([x, y, *z])
sos = 10 * sum((zi + (x + y) / 5) ** 2 for zi in z) + sp.Rational(7, 2) * (x + y) ** 2 + sp.Rational(9, 2) * (x - y) ** 2
if sp.expand((v.T * P * v)[0] - sos) != 0: fails.append("SOS identity")
D = sp.diag(1, 20, 2, 2, 2, 2, 2); M = P * D * P.inv(); d = [M[i, i] for i in range(7)]
print("diag(PDP^-1) =", d, "sum", sum(d), "| spectrum (1,20,2,2,2,2,2)")
print("min diag", min(d), "> min eigenvalue 1:", min(d) > 1, "| max diag", max(d), "> 20:", max(d) > 20)
if not (min(d) > 1 and max(d) > 20 and sum(d) == 31): fails.append("witness")
c, r = sp.symbols('c r', positive=True)
Mcr = P * sp.diag(c, c + r, c + 1, c + 1, c + 1, c + 1, c + 1) * P.inv()
fam = [sp.simplify(Mcr[i, i] - c) for i in range(7)]
print("family diag - c:", fam)
if fam != [(r - 18) / 63, (80 * r - 18) / 63] + [(39 - 2 * r) / 35] * 5: fails.append("family")
# control: n <= 4 (proven range in the source): random PD examples satisfy min diag <= min eig <= max eig <= max diag
import random; random.seed(1)
for _ in range(200):
    A = sp.Matrix(4, 4, lambda i, j: random.randint(-3, 3)); Pc = A * A.T + sp.eye(4)
    dd = [random.randint(1, 9) for _ in range(4)]; Mc = Pc * sp.diag(*dd) * Pc.inv(); dg = [Mc[i, i] for i in range(4)]
    if not (min(dg) <= min(dd) and max(dg) >= max(dd)): fails.append("control-n4"); break
print("control n=4: necessary majorization inequalities hold on 200 random PD examples:", "control-n4" not in fails)
print("RESULT", "PASS" if not fails else f"FAIL {fails}")
sys.exit(0 if not fails else 1)
