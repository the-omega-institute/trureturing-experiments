"""Numerical check of a counterfamily to the open question of Johnston, Moein, Pereira and Plosker,
"Absolutely k-incoherent quantum states and spectral inequalities for factor width of a matrix",
arXiv:2205.05110v2, Phys. Rev. A 106, 052417 (2022), Section 5: is the purity condition of Theorem 7,
Tr(rho^2) <= 1/(n-1), necessary for absolute 2-incoherence when n >= 4?

Proposed answer (n >= 5): no.  With r = sqrt((n-1)(n-2)/2), a = (1+r)/n, b = (1 - 2r/(n-2))/n,
B = diag(a, a, b, ..., b), rho0 = (I - B)/(n-1) and rho_t = (1+t) rho0 - (t/n) I, the state rho_t is
absolutely 2-incoherent for small t > 0 while Tr(rho_t^2) = 1/n + (1+t)^2/(n(n-1)) > 1/(n-1).
The proof rests on a dual gap: 1 - <A, U B U*> >= h0 > 0 for every A with all 2x2 principal minors PSD
(the cone L2) and Tr A = 1, and every unitary U.

Prediction registered before the run:
  (1) the spectral identities 2a + (n-2)b = 1, 2a^2 + (n-2)b^2 = 1 hold and rho_t has the stated purity;
  (2) for n = 5, 6, the largest <A, U B U*> found (random unitaries plus local search) stays below
      1 - 0.005;
  (3) for n = 5 and t = 1e-3 (far above the proof's t), U rho_t U* has an explicit factor-width-2
      decomposition (an SDP with 2x2 blocks) for every tested unitary, including the adversarial ones.
Negative controls: (a) at n = 4 the same construction has gap below 1e-3 (the orbit reaches L2, which is
why n = 4 is excluded); (b) a pure state is not absolutely 2-incoherent (the decomposition SDP is
infeasible for some unitary); (c) Theorem 7's sufficient side: a state with Tr(rho^2) < 1/(n-1) decomposes.
Usage: python3 check.py  (needs numpy and cvxpy with the Clarabel solver; seed fixed)
"""
import sys, itertools
import numpy as np
import cvxpy as cp

rng = np.random.default_rng(20261017)
OK = True


def fail(m):
    global OK
    OK = False
    print("FAIL", m)


def construction(n):
    r = np.sqrt((n - 1) * (n - 2) / 2)
    a = (1 + r) / n
    b = (1 - 2 * r / (n - 2)) / n
    return a, b, np.diag([a, a] + [b] * (n - 2))


def haar(n):
    Z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / abs(np.diag(R)))


def inner_max(Bp):
    n = Bp.shape[0]
    A = cp.Variable((n, n), hermitian=True)
    cons = [cp.real(cp.trace(A)) == 1]
    for i, j in itertools.combinations(range(n), 2):
        cons.append(cp.bmat([[A[i, i], A[i, j]], [A[j, i], A[j, j]]]) >> 0)
    pr = cp.Problem(cp.Maximize(cp.real(cp.trace(A @ Bp))), cons)
    pr.solve(solver=cp.CLARABEL)
    return pr.value


def fw2_slack(X):
    """largest s such that X - s*I/n is a sum of PSD blocks supported on coordinate pairs (s >= 0 means X in FW2)."""
    n = X.shape[0]
    blocks = {}
    expr = 0
    for i, j in itertools.combinations(range(n), 2):
        P = cp.Variable((2, 2), hermitian=True)
        blocks[(i, j)] = P
        E = np.zeros((n, 2))
        E[i, 0] = 1
        E[j, 1] = 1
        expr = expr + E @ P @ E.T
    s = cp.Variable()
    cons = [P >> 0 for P in blocks.values()] + [expr == X - s * np.eye(n) / n]
    pr = cp.Problem(cp.Maximize(s), cons)
    pr.solve(solver=cp.CLARABEL)
    return pr.value


def max_gap_search(B, restarts=4, iters=250):
    n = B.shape[0]
    best = -1
    for _ in range(40):
        U = haar(n)
        best = max(best, inner_max(U @ B @ U.conj().T))
    worstU = None
    for _ in range(restarts):
        U = haar(n)
        cur = inner_max(U @ B @ U.conj().T)
        for it in range(iters):
            H = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            H = (H + H.conj().T) / 2
            w, V = np.linalg.eigh(H)
            step = 0.3 * np.exp(-it / 80)
            Ut = (V @ np.diag(np.exp(1j * step * w)) @ V.conj().T) @ U
            v = inner_max(Ut @ B @ Ut.conj().T)
            if v > cur:
                cur, U = v, Ut
        if cur > best:
            best, worstU = cur, U
    return best, worstU


print("== (1) spectral identities and purity")
for n in range(5, 9):
    a, b, B = construction(n)
    if not (abs(2 * a + (n - 2) * b - 1) < 1e-12 and abs(2 * a * a + (n - 2) * b * b - 1) < 1e-12 and a > 0 > b):
        fail(f"identities n={n}")
    t = 1 / (100000 * n ** 4)
    rho = (1 + t) * (np.eye(n) - B) / (n - 1) - t / n * np.eye(n)
    pur = np.trace(rho @ rho).real
    form = 1 / n + (1 + t) ** 2 / (n * (n - 1))
    print(f"  n={n}: a={a:.6f} b={b:.6f} min eig rho_t={np.linalg.eigvalsh(rho).min():.6f} purity-1/(n-1)={pur - 1/(n-1):.3e} formula match {abs(pur-form)<1e-14}")
    if not (np.linalg.eigvalsh(rho).min() > 0 and pur > 1 / (n - 1) and abs(pur - form) < 1e-14):
        fail(f"purity n={n}")

print("== (2) dual gap search: max <A, U B U*> over A in L2, Tr A = 1")
worst = {}
for n in (5, 6):
    _, _, B = construction(n)
    m, U = max_gap_search(B)
    worst[n] = U
    print(f"  n={n}: largest value found {m:.6f}, gap {1-m:.6f}")
    if not m < 1 - 0.005:
        fail(f"gap n={n}")

print("== (3) primal decompositions at n = 5, t = 1e-3")
n = 5
_, _, B = construction(n)
t = 1e-3
rho = (1 + t) * (np.eye(n) - B) / (n - 1) - t / n * np.eye(n)
print(f"  purity {np.trace(rho @ rho).real:.6f} > 1/(n-1) = {1/(n-1):.6f}")
tests = [haar(n) for _ in range(30)] + ([worst[5]] if worst.get(5) is not None else [])
slacks = [fw2_slack(U @ rho @ U.conj().T) for U in tests]
print(f"  min slack over {len(tests)} unitaries (incl. the adversarial one): {min(slacks):.3e}")
if not min(slacks) > 0:
    fail("primal decomposition")

print("== negative controls")
a4, b4, B4 = construction(4)
m4, _ = max_gap_search(B4, restarts=3, iters=200)
print(f"  (a) n=4: largest <A, U B U*> {m4:.6f}, gap {1-m4:.2e}")
if not 1 - m4 < 1e-3:
    fail("n=4 control")
psi = np.zeros(5, complex)
psi[0] = 1
pure = np.outer(psi, psi.conj())
sl = min(fw2_slack(U @ pure @ U.conj().T) for U in [haar(5) for _ in range(5)])
print(f"  (b) pure state: min slack {sl:.3e} (negative: not 2-incoherent after a generic unitary)")
if not sl < -1e-3:
    fail("pure control")
mixed = 0.9 * np.eye(5) / 5 + 0.1 * pure
print(f"  (c) purity {np.trace(mixed @ mixed).real:.4f} < 1/4: min slack {min(fw2_slack(U @ mixed @ U.conj().T) for U in [haar(5) for _ in range(10)]):.3e}")
print("ALL_OK" if OK else "SOME_FAILED")
sys.exit(0 if OK else 1)
