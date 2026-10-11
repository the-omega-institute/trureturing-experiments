"""cHaar moment-operator norm under concatenation.

Source: Duschenes, García-Martín, Holmes, Cerezo, "Moments of Quantum Channel Ensembles",
arXiv:2511.12700v2 (PRX Quantum 7, 033050 (2026)), Appendix B, Conjecture 1, second clause:
||T^(t)k'||_HS^2 <= ||T^(t)k||_HS^2 for all k' >= k.

Checks, for t = 2:
 1. The second moment operator T_2 of the channel Lambda_U(X) = Tr_E[U (X (x) |0><0|) U^dag], U Haar on
    U(de), built EXACTLY (sympy rationals) from the two-copy Haar twirl formula
    P(Y) = [Tr Y - Tr(F Y)/D] I/(D^2-1) + [Tr(F Y) - Tr Y/D] F/(D^2-1), D = de,
    as a d^4 x d^4 Liouville matrix acting on X in B(H (x) H).
 2. N_k = ||T_2^k||_HS^2 from exact matrix powers equals the closed form
    1 + q (1 - lam^k)^2 + lam^(2k), lam = e(d^2-1)/(d^2 e^2 - 1), q = (d^2-1)/(d^2 e + 1)^2,
    for 2 <= d, e <= 3 and 1 <= k <= 8.
 3. N_5 > N_4 at d = e = 2 (exact): N_4 = 1213852/1171875, N_5 = 30362908/29296875.
 4. For every 2 <= d, e <= 6 the sequence N_k increases strictly from some k on (exact closed form, k <= 60).
 5. Monte Carlo (seeded, 20000 Haar samples) of T_2 at d = e = 2 agrees with the exact matrix to 0.02.
"""
import itertools, sys
import numpy as np
import sympy as sp

bad = 0
def ok(c, msg):
    global bad
    print(("ok   " if c else "FAIL ") + msg)
    if not c: bad += 1

def T2_exact(d, e):
    D = d * e
    idx = list(itertools.product(range(d), repeat=2))  # basis of H(x)H, pairs (a1,a2)
    n = len(idx)
    # Liouville matrix L[(r,c),(r',c')] = T2(E_{r' c'})[r,c], with r,c in idx
    L = sp.zeros(n * n, n * n)
    for ri, ci in itertools.product(range(n), repeat=2):
        X = sp.zeros(n, n); X[ri, ci] = 1
        # Y = X (x) nu (x) nu reordered; Tr(Y) = Tr X; Tr(F_D Y) = Tr(F X) (purity of nu)
        trX = X.trace()
        F = sp.zeros(n, n)
        for i, (a, b) in enumerate(idx):
            F[idx.index((b, a)), i] = 1
        trFX = (F * X).trace()
        # partial traces over environment: Tr_E I_D(x)I_D = e^2 I, Tr_E F_D = e F
        Id = sp.eye(n)
        out = (sp.Rational(e * e) * (trX - trFX / D) / (D * D - 1)) * Id + (sp.Rational(e) * (trFX - trX / D) / (D * D - 1)) * F
        for (oi, o2) in itertools.product(range(n), repeat=2):
            L[oi * n + o2, ri * n + ci] = out[oi, o2]
    return L

def Nclosed(d, e, k):
    lam = sp.Rational(e * (d * d - 1), d * d * e * e - 1)
    q = sp.Rational(d * d - 1, (d * d * e + 1) ** 2)
    return 1 + q * (1 - lam ** k) ** 2 + lam ** (2 * k)

for d, e in [(2, 2), (2, 3), (3, 2)]:
    L = T2_exact(d, e)
    P = sp.eye(L.shape[0])
    allok = True
    for k in range(1, 9):
        P = P * L
        Nk = sum(x * x for x in P)  # real rational entries
        if sp.simplify(Nk - Nclosed(d, e, k)) != 0:
            allok = False
    ok(allok, f"exact matrix powers match the closed form, d={d}, e={e}, k<=8")

N4, N5 = Nclosed(2, 2, 4), Nclosed(2, 2, 5)
ok(N4 == sp.Rational(1213852, 1171875) and N5 == sp.Rational(30362908, 29296875), "N_4, N_5 at d=e=2")
ok(N5 - N4 == sp.Rational(5536, 9765625) and N5 > N4, "N_5 - N_4 = 5536/9765625 > 0")
for d in range(2, 7):
    for e in range(2, 7):
        seq = [Nclosed(d, e, k) for k in range(1, 61)]
        k0 = next((k for k in range(len(seq) - 1) if all(seq[j + 1] > seq[j] for j in range(k, len(seq) - 1))), None)
        ok(k0 is not None, f"d={d}, e={e}: strictly increasing from k={None if k0 is None else k0+1} up to k=60")

# Monte Carlo
rng = np.random.default_rng(20261012)
d, e = 2, 2; D = d * e
def haar(n):
    z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    q, r = np.linalg.qr(z); return q * (np.diag(r) / abs(np.diag(r)))
def liouville_two_copy(U):
    U4 = U.reshape(d, e, d, e)
    K = [U4[:, j, :, 0] for j in range(e)]
    S = np.zeros((d ** 4, d ** 4), dtype=complex)
    for A in K:
        for B in K:
            M = np.kron(A, B); S += np.kron(M, M.conj())
    return S
M = 20000
T = sum(liouville_two_copy(haar(D)) for _ in range(M)) / M
Lx = np.array(T2_exact(2, 2)).astype(float)
ok(np.max(abs(T - Lx)) < 0.02, f"Monte Carlo agrees with the exact T_2 (max deviation {np.max(abs(T - Lx)):.4f})")
print("bad=", bad)
sys.exit(1 if bad else 0)
