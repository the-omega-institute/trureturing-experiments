"""Checks for Saleem, Shaji, Babu, Luo, Langfitt, Yu and Gray, arXiv:2504.13729v4, Section IV.B: "in all cases we
conjecture that CoE <= QFI with equality at specific points" for the open evolution (41) (flip-flop Hamiltonian
H_g = g(|01><10| + |10><01|), amplitude damping sigma_- = |1><0| on both qubits at rate kappa, Appendix B convention).
(1) Exact (sympy) family (the X-block spin-flip spectrum is checked symbolically for every block state): rho_0 = diag(1, 11, 3, 1)/16; the explicit rho_g(t) solves (41) for every g, kappa, t;
    it is full rank; at g* = pi/(4t) the concurrence is on a smooth positive branch with CoE = 2 q t^2 and the SLD
    gives F_Q = 8 q t^2/(8 - q), q = exp(-kappa t); CoE - F_Q = 2 q t^2 (4 - q)/(8 - q) > 0.
    Explicit point t = 1, kappa = ln 2, g = pi/4: C = (8 - sqrt 33)/32, F_Q = 8/15, CoE = 1.
(2) Independent numerics (numpy): the Liouvillian exponential reproduces (1) and the source's Eqs. (42)-(43) for
    rho_0 = |01><01|; with the same code, kappa = 0 shows no
    violation (the source's pure-state theorem), while the pure product initial state (|00> + |10>)/sqrt 2 violates
    the inequality at smooth points (e.g. kappa = 0.3, t = 3, g = 2.9).  Prints ALL_OK."""
import sympy as sp, numpy as np
g, t, k = sp.symbols('g t kappa', positive=True)
q = sp.exp(-k*t); th = 2*g*t
A = 1 - q + q**2/16; D = q**2/16; K = 1 - q/8; Dl = sp.Rational(1, 2)
rho = sp.Matrix([[D, 0, 0, 0], [0, q*(K + Dl*sp.cos(th))/2, sp.I*q*Dl*sp.sin(th)/2, 0],
                 [0, -sp.I*q*Dl*sp.sin(th)/2, q*(K - Dl*sp.cos(th))/2, 0], [0, 0, 0, A]])
# basis order |00>, |01>, |10>, |11>; sigma_- = |1><0| on each qubit
sm = sp.Matrix([[0, 0], [1, 0]]); I2 = sp.eye(2)
def kron(a, b): return sp.Matrix(sp.BlockMatrix([[a[i, j]*b for j in range(a.cols)] for i in range(a.rows)]).as_explicit())
J1, J2 = kron(sm, I2), kron(I2, sm)
X = sp.zeros(4); X[1, 2] = X[2, 1] = 1
H = g*X
def Dis(J, r): return J*r*J.H - (J.H*J*r + r*J.H*J)/2
L = -sp.I*(H*rho - rho*H) + k*(Dis(J1, rho) + Dis(J2, rho))
ok1 = sp.simplify(sp.diff(rho, t) - L) == sp.zeros(4)
ok1 &= sp.simplify(rho.subs(t, 0) - sp.diag(1, 11, 3, 1)/16) == sp.zeros(4) and sp.simplify(rho.trace() - 1) == 0
# explicit point
pt = {t: 1, k: sp.log(2), g: sp.pi/4}
r1 = sp.nsimplify(sp.simplify(rho.subs(pt)))
ok1 &= r1 == sp.Matrix([[1, 0, 0, 0], [0, 15, 8*sp.I, 0], [0, -8*sp.I, 15, 0], [0, 0, 0, 33]])/64
YY = kron(sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]))
R = r1*YY*r1.conjugate()*YY
ev = sorted([sp.nsimplify(e) for e in R.eigenvals(multiple=True)], reverse=True)
ok1 &= ev == [sp.Rational(529, 4096), sp.Rational(49, 4096), sp.Rational(33, 4096), sp.Rational(33, 4096)]
lam = [sp.sqrt(e) for e in ev]
Cpt = sp.simplify(lam[0] - lam[1] - lam[2] - lam[3]); ok1 &= sp.simplify(Cpt - (8 - sp.sqrt(33))/32) == 0 and Cpt > 0
# general block form: the spin-flip product has characteristic polynomial (x - AD)^2 (x^2 - 2(bc + y^2)x + (bc - y^2)^2),
# i.e. square roots sqrt(AD), sqrt(AD), sqrt(bc) + |y|, sqrt(bc) - |y|, so C = 2(|y| - sqrt(AD)) whenever |y| > sqrt(AD)
aS, dS, bS, cS = sp.symbols('A_ D_ b_ c_', positive=True); yS = sp.Symbol('y_', real=True); xS = sp.Symbol('x_')
rX = sp.Matrix([[dS, 0, 0, 0], [0, bS, sp.I*yS, 0], [0, -sp.I*yS, cS, 0], [0, 0, 0, aS]])
cpX = (rX*YY*rX.conjugate()*YY).charpoly(xS).as_expr()
ok1 &= sp.expand(cpX - (xS - aS*dS)**2*(xS**2 - 2*(bS*cS + yS**2)*xS + (bS*cS - yS**2)**2)) == 0
ok1 &= sp.expand(((sp.sqrt(bS*cS) + yS)**2)*((sp.sqrt(bS*cS) - yS)**2) - (bS*cS - yS**2)**2) == 0
# concurrence branch near g*: C(g) = 2(|z| - sqrt(A D)) with z = q sin(th)/4 > 0 near g*
Cg = 2*(q*sp.sin(th)/4 - sp.sqrt(A*D))
CoE = sp.simplify(-sp.diff(Cg, g, 2).subs(g, sp.pi/(4*t)))
ok1 &= sp.simplify(CoE - 2*q*t**2) == 0
# SLD on the full-rank state at g*: L = 0 on |00>,|11>; (2 t Dl / K) V on the middle block
V = sp.Matrix([[-sp.sin(th), sp.I*sp.cos(th)], [-sp.I*sp.cos(th), sp.sin(th)]])  # -sin(th) Z - cos(th) Y
Lsld = sp.zeros(4); Lsld[1:3, 1:3] = (2*t*Dl/K)*V
dr = sp.diff(rho, g)
ok1 &= sp.simplify((rho*Lsld + Lsld*rho)/2 - dr) == sp.zeros(4)
FQ = sp.simplify((rho*Lsld*Lsld).trace())
ok1 &= sp.simplify(FQ - 8*q*t**2/(8 - q)) == 0
gapf = sp.simplify(CoE - FQ - 2*q*t**2*(4 - q)/(8 - q)); ok1 &= gapf == 0
ok1 &= sp.nsimplify(FQ.subs(pt)) == sp.Rational(8, 15) and sp.nsimplify(CoE.subs(pt)) == 1
# the branch conditions at g*: sin(th) = 1 and q/4 > q sqrt(A)/4 (i.e. A < 1) for 0 < q < 1
print("exact family identities (ODE, trace, point state, R spectrum, concurrence, SLD, F_Q, CoE - F_Q):", ok1)
# (2) numerics
def nexpm(Am):
    n = np.linalg.norm(Am, 1); s = max(0, int(np.ceil(np.log2(n))) + 1) if n > 0 else 0
    B = Am/2**s; E = np.eye(Am.shape[0], dtype=complex); T = E.copy()
    for j in range(1, 25): T = T@B/j; E = E + T
    for _ in range(s): E = E@E
    return E
smn = np.array([[0, 0], [1, 0]], dtype=complex); Jn = [np.kron(smn, np.eye(2)), np.kron(np.eye(2), smn)]
Xn = np.zeros((4, 4)); Xn[1, 2] = Xn[2, 1] = 1
Yn = np.array([[0, -1j], [1j, 0]]); YYn = np.kron(Yn, Yn)
def evolve(r0, gg, kk, tt):
    I4 = np.eye(4); Hn = gg*Xn
    Lv = -1j*(np.kron(I4, Hn) - np.kron(Hn.T, I4))
    for Jm in Jn:
        JdJ = Jm.conj().T@Jm
        Lv = Lv + kk*(np.kron(Jm.conj(), Jm) - 0.5*np.kron(I4, JdJ) - 0.5*np.kron(JdJ.T, I4))
    return (nexpm(Lv*tt)@r0.reshape(-1, order='F')).reshape(4, 4, order='F')
def conc(r):
    e = np.sort(np.sqrt(np.abs(np.linalg.eigvals(r@YYn@r.conj()@YYn))))[::-1]
    return max(0.0, e[0] - e[1] - e[2] - e[3]), e
def qfi(r, d):
    w, U = np.linalg.eigh(r); dd = U.conj().T@d@U
    return sum(2*abs(dd[i, j])**2/(w[i] + w[j]) for i in range(4) for j in range(4) if w[i] + w[j] > 1e-12)
def gap(r0, gg, kk, tt, h):
    rp, r0g, rm = [evolve(r0, gg + e, kk, tt) for e in (h, 0, -h)]
    coe = -(conc(rp)[0] - 2*conc(r0g)[0] + conc(rm)[0])/h**2
    return coe - qfi(r0g, (rp - rm)/(2*h)), conc(r0g)
r0 = np.diag([1, 11, 3, 1]).astype(complex)/16
num = [gap(r0, np.pi/4, np.log(2), 1.0, h)[0] for h in (1e-3, 1e-4)]
ok2 = all(abs(x - 7/15) < 1e-5 for x in num)
rng = np.random.default_rng(0); worst = -1e9; cnt = 0
for _ in range(3000):
    psi = rng.normal(size=4) + 1j*rng.normal(size=4); psi /= np.linalg.norm(psi)
    rr = np.outer(psi, psi.conj()); tt = rng.uniform(0.1, 4); gg = rng.uniform(0.05, 3)
    c, e = conc(evolve(rr, gg, 0.0, tt))
    if c < 0.05: continue
    ds = [gap(rr, gg, 0.0, tt, h)[0] for h in (1e-3, 3e-4, 1e-4)]
    if max(ds) - min(ds) > 1e-3*max(1, abs(ds[0])): continue
    cnt += 1; worst = max(worst, ds[-1])
ok2 &= cnt > 0 and worst < 1e-6
psi = np.array([1, 0, 1, 0], dtype=complex)/np.sqrt(2); rr = np.outer(psi, psi.conj())
ds = [gap(rr, 2.9, 0.3, 3.0, h) for h in (1e-3, 3e-4, 1e-4)]
c, e = ds[-1][1]
ok3 = all(d[0] > 3 for d in ds) and max(d[0] for d in ds) - min(d[0] for d in ds) < 1e-3 and c > 0.04 and e[0] - e[1] > 0.04
# the conventions reproduce the source's Eqs. (42)-(43) for rho_opt(0) = |01><01|: F = 4 t^2 e^{-kappa t}, C = e^{-kappa t}|sin 2gt|
r01 = np.zeros((4, 4), dtype=complex); r01[1, 1] = 1
src_ok = True
for (gg, kk, tt) in [(0.4, 0.5, 1.3), (1.1, 0.2, 2.0), (0.7, 1.0, 0.6)]:
    rg = evolve(r01, gg, kk, tt); h = 1e-5
    Fq = qfi(rg, (evolve(r01, gg + h, kk, tt) - evolve(r01, gg - h, kk, tt))/(2*h))
    src_ok &= abs(Fq - 4*tt**2*np.exp(-kk*tt)) < 1e-6 and abs(conc(rg)[0] - np.exp(-kk*tt)*abs(np.sin(2*gg*tt))) < 1e-7
ok2 &= src_ok
print("numerics: family point gap", [round(x, 6) for x in num], "; source Eqs. (42)-(43) reproduced", src_ok, "; kappa = 0 smooth samples", cnt, "max gap", round(worst, 4),
      "; pure product state (|00>+|10>)/sqrt2 at kappa=0.3, t=3, g=2.9 gap", round(ds[-1][0], 4), "C", round(c, 4))
print("ALL_OK" if ok1 and ok2 and ok3 else "FAIL")
