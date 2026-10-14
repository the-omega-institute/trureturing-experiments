#!/usr/bin/env python3
"""Checks for Campbell, "Optimal entangling capacity of dynamical processes", arXiv:1007.1445v2 (PRA 82, 042314 (2010)),
Section VII: "whether log-negativity EC is resource independent for all unitaries, or indeed all CPTP maps".
Candidate separation (all-CPTP clause only): Q = I_9 - P, P the projector onto the five Tiles product vectors
(Bennett et al., PRL 82, 5385 (1999)); Lambda_t(X) = t Tr(QX) Phi + Tr((I - tQ)X) |22><22| on two qutrits.
1. Exact (sympy): the Tiles vectors are orthonormal; Q^2 = Q, rank Q = 4, Q^{Gamma_B} = Q.
2. Exact (sympy): the ten 3x3 certificates (determinants, squared inverse Frobenius norms) and hence eps >= 1/6561.
3. Numerics: beta = max over product states of <xy|Q|xy> (alternating maximisation, many restarts) satisfies beta <= 1 - 1/6561.
4. Numerics: the Kraus family {sqrt(t)|phi><f_j|Q, |22><f_j|(P + sqrt(1-t)Q)} is complete and implements Lambda_t;
   ||(Lambda_t(Q/4))^{Gamma}||_1 = 1 + t and Q/4 is PPT (log-negativity 0).
5. Numerics: random separable inputs with qubit ancillas on both sides give ||Gamma(Lambda_t (x) id)(sigma)||_1 <= 1 + t*beta_num.
Prints ALL_OK.
"""
import itertools, numpy as np, sympy as sp
ok = True
e = [sp.Matrix([1,0,0]), sp.Matrix([0,1,0]), sp.Matrix([0,0,1])]
a = [e[0], e[0]-e[1], e[2], e[1]-e[2], e[0]+e[1]+e[2]]
b = [e[0]-e[1], e[2], e[1]-e[2], e[0], e[0]+e[1]+e[2]]
def kron(u, v): return sp.Matrix([u[i]*v[j] for i in range(3) for j in range(3)])
u = [kron(a[i], b[i]) / sp.sqrt((a[i].T*a[i])[0] * (b[i].T*b[i])[0]) for i in range(5)]
G = sp.Matrix(5, 5, lambda i, j: sp.simplify((u[i].T*u[j])[0]))
if G != sp.eye(5): print('GRAM_FAIL'); ok = False
P = sum((ui*ui.T for ui in u), sp.zeros(9, 9)); Q = sp.eye(9) - P
if sp.simplify(Q*Q - Q) != sp.zeros(9, 9) or Q.rank() != 4: print('PROJ_FAIL'); ok = False
def ptB(M, dA=3, dB=3):
    R = sp.zeros(dA*dB, dA*dB)
    for i, j, k, l in itertools.product(range(dA), range(dB), range(dA), range(dB)):
        R[i*dB+j, k*dB+l] = M[i*dB+l, k*dB+j]
    return R
if sp.simplify(ptB(Q) - Q) != sp.zeros(9, 9): print('PT_FAIL'); ok = False
print('exact: Tiles orthonormal, Q^2=Q, rank 4, Q^Gamma=Q')
dets, fro = [], []
for I in itertools.combinations(range(5), 3):
    C = sp.Matrix([list(a[i].T) for i in I])
    dets.append(C.det()); Ci = C.inv(); fro.append(sum(x**2 for x in Ci))
print('triples', list(itertools.combinations(range(5), 3)))
print('dets', dets, 'squared inverse Frobenius', fro)
if any(f > 9 for f in fro): print('CERT_FAIL'); ok = False
# eps >= 1/6561: each triple gives sum |<a_i,x>|^2 >= 1/||C^-1||_F^2 >= 1/9 -> three overlaps >= 1/27; intersecting index; ||a||^2||b||^2 <= 9
nrm = [(a[i].T*a[i])[0]*(b[i].T*b[i])[0] for i in range(5)]
if max(nrm) > 9: print('NORM_FAIL'); ok = False
print('exact certificate: eps >= (1/27)^2/9 =', sp.Rational(1, 27)**2 / 9)
Qn = np.array(Q.evalf(), dtype=float); Pn = np.eye(9) - Qn
rng = np.random.default_rng(1)
best = 0.0
for _ in range(3000):
    x = rng.normal(size=3) + 1j*rng.normal(size=3); x /= np.linalg.norm(x)
    for _ in range(60):
        My = np.einsum('i,ijkl,k->jl', x.conj(), Qn.reshape(3,3,3,3), x); w, V = np.linalg.eigh(My); y = V[:, -1]
        Mx = np.einsum('j,ijkl,l->ik', y.conj(), Qn.reshape(3,3,3,3), y); w, V = np.linalg.eigh(Mx); x = V[:, -1]
    best = max(best, float(np.real(np.kron(x, y).conj() @ Qn @ np.kron(x, y))))
print('beta_num = %.10f   1 - beta_num = %.3e   certificate bound 1 - 1/6561 = %.10f' % (best, 1-best, 1-1/6561))
if best > 1 - 1/6561 + 1e-12: print('BETA_FAIL'); ok = False
phi = np.zeros(9); phi[0] = phi[4] = 1/np.sqrt(2); Phi = np.outer(phi, phi); F = np.zeros((9, 9)); F[8, 8] = 1
def tn(M): return float(np.sum(np.abs(np.linalg.eigvalsh((M + M.conj().T)/2))))
def ptBn(M, dA, dB):
    return M.reshape(dA, dB, dA, dB).transpose(0, 3, 2, 1).reshape(dA*dB, dA*dB)
for t in (0.25, 0.5, 0.9):
    K = [np.sqrt(t)*np.outer(phi, np.eye(9)[j]) @ Qn for j in range(9)] + \
        [np.outer(np.eye(9)[8], np.eye(9)[j]) @ (Pn + np.sqrt(1-t)*Qn) for j in range(9)]
    comp = sum(k.conj().T @ k for k in K)
    if np.abs(comp - np.eye(9)).max() > 1e-12: print('KRAUS_FAIL'); ok = False
    X = rng.normal(size=(9, 9)) + 1j*rng.normal(size=(9, 9)); X = X @ X.conj().T; X /= np.trace(X)
    lam = t*np.trace(Qn @ X)*Phi + np.trace((np.eye(9) - t*Qn) @ X)*F
    if np.abs(sum(k @ X @ k.conj().T for k in K) - lam).max() > 1e-12: print('ACTION_FAIL'); ok = False
    tau = Qn/4; out = t*Phi + (1-t)*F
    if abs(tn(ptBn(tau, 3, 3)) - 1) > 1e-12 or abs(tn(ptBn(out, 3, 3)) - (1+t)) > 1e-12: print('ASSIST_FAIL'); ok = False
    # separable inputs with qubit ancillas: sigma on (A A') (B B'), A=B=3, A'=B'=2; apply Lambda on (A,B)
    worst = 0.0
    for _ in range(400):
        sig = np.zeros((36, 36), dtype=complex); ps = rng.dirichlet(np.ones(4))
        for p in ps:
            xa = rng.normal(size=6) + 1j*rng.normal(size=6); xa /= np.linalg.norm(xa)
            yb = rng.normal(size=6) + 1j*rng.normal(size=6); yb /= np.linalg.norm(yb)
            v = np.kron(xa, yb); sig += p*np.outer(v, v.conj())
        S = sig.reshape(3, 2, 3, 2, 3, 2, 3, 2)               # (a a' b b') x (a a' b b')
        Sab = S.transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(9, 4, 9, 4)  # (a b)(a' b') x (a b)(a' b')
        outk = np.zeros((9, 4, 9, 4), dtype=complex)
        for k in K: outk += np.einsum('ij,jakb,lk->ialb', k, Sab, k.conj())
        O = outk.reshape(3, 3, 2, 2, 3, 3, 2, 2).transpose(0, 2, 1, 3, 4, 6, 5, 7).reshape(36, 36)
        worst = max(worst, tn(ptBn(O, 6, 6)))
    print('t=%.2f  assisted norm %.6f (=1+t)  max separable norm with qubit ancillas %.6f  bound 1+t*beta %.6f' % (t, tn(ptBn(out, 3, 3)), worst, 1+t*best))
    if worst > 1 + t*best + 1e-9: print('SEP_FAIL'); ok = False
print('ALL_OK' if ok else 'FAIL')
