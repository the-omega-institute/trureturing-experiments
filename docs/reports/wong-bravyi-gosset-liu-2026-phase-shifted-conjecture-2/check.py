"""Wong-Bravyi-Gosset-Liu, arXiv:2602.03605v1, Section 5.6: does Conjecture 2 (the second smallest eigenvalue
of H_s equals the smallest eigenvalue of H_s on the odd-weight subspace) hold for phase-shifted EPR-like
Hamiltonians H_s(theta) = - sum_{(p,q) in E} h_{p,q}(s e^{i theta_{p,q}})? Counterexample: the 4-cycle with phases
(0, 0, 0, pi), s = 3/5 (exact, with sympy), plus a scan."""
import sys, itertools
import numpy as np
import sympy as sp

ok = True
def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name); ok = ok and bool(cond)

def ham(n, edges, s, thetas, exact=False):
    N = 2 ** n
    H = sp.zeros(N, N) if exact else np.zeros((N, N), dtype=complex)
    for (p, q), th in zip(edges, thetas):
        z = s * (sp.exp(sp.I * th) if exact else np.exp(1j * th))
        zc = sp.conjugate(z) if exact else np.conj(z)
        for x in range(N):
            bp, bq = (x >> (n - 1 - p)) & 1, (x >> (n - 1 - q)) & 1
            flip = x ^ (1 << (n - 1 - p)) ^ (1 << (n - 1 - q))
            if bp == 0 and bq == 0:
                H[x, x] -= 1; H[flip, x] -= z          # z |11><00|
            elif bp == 1 and bq == 1:
                H[x, x] -= z * zc; H[flip, x] -= zc    # z* |00><11|
    return H

edges = [(0, 1), (1, 2), (2, 3), (3, 0)]
s = sp.Rational(3, 5)
H = ham(4, edges, s, [0, 0, 0, sp.pi], exact=True)
check("H is Hermitian (exact)", sp.simplify(H - H.H) == sp.zeros(16, 16))
odd = [x for x in range(16) if bin(x).count("1") % 2 == 1]; even = [x for x in range(16) if bin(x).count("1") % 2 == 0]
check("H preserves Hamming-weight parity", all(H[i, j] == 0 for i in odd for j in even))
Ho = H.extract(odd, odd); He = H.extract(even, even)
po = sp.factor(Ho.charpoly().as_expr()); pe = sp.factor(He.charpoly().as_expr())
print("odd charpoly:", po); print("even charpoly:", pe)
ev = sorted([complex(v).real for v in np.linalg.eigvalsh(np.array(H.evalf(30).tolist(), dtype=complex))])
evo = sorted(np.linalg.eigvalsh(np.array(Ho.evalf(30).tolist(), dtype=complex)).real)
print("spectrum:", np.round(ev, 6)); print("odd minimum:", round(evo[0], 6))
# exact comparison: lambda_2 = -beta_+ with t = s^2, odd minimum = -omega_+
t = s ** 2
beta_p = (1 + 5 * t + sp.sqrt((1 + t) * (1 + 9 * t))) / 2
omega_p = 1 + t + sp.sqrt(1 + t ** 2)
alpha_p = (5 + t + sp.sqrt((1 + t) * (9 + t))) / 2
lam = sp.Symbol('x')
ce = He.charpoly(lam).as_expr(); co = Ho.charpoly(lam).as_expr()
check("exact: -alpha_+ and -beta_+ are roots of the even characteristic polynomial and -omega_+ of the odd one",
      all(sp.simplify(sp.expand(e)) == 0 for e in (ce.subs(lam, -alpha_p), ce.subs(lam, -beta_p), co.subs(lam, -omega_p))))
check(f"exact: beta_+ - omega_+ = {sp.nsimplify(beta_p - omega_p)} > 0 at s = 3/5", sp.N(beta_p - omega_p, 50) > 0)
check(f"numeric: second smallest eigenvalue {ev[1]:.6f} < odd minimum {evo[0]:.6f}, lambda_1 < lambda_2 simple", ev[1] < evo[0] - 1e-9 and ev[0] < ev[1] - 1e-9)
check("numeric: lambda_1 = -alpha_+, lambda_2 = -beta_+, odd minimum = -omega_+",
      abs(ev[0] + float(alpha_p)) < 1e-9 and abs(ev[1] + float(beta_p)) < 1e-9 and abs(evo[0] + float(omega_p)) < 1e-9)
# unphased C4 at the same s satisfies Conjecture 2
H0 = ham(4, edges, 0.6, [0, 0, 0, 0]); e0 = sorted(np.linalg.eigvalsh(H0)); o0 = sorted(np.linalg.eigvalsh(H0[np.ix_(odd, odd)]))
check(f"unphased C4 at s = 0.6: lambda_2 {e0[1]:.6f} equals odd minimum {o0[0]:.6f}", abs(e0[1] - o0[0]) < 1e-9)
# scan of s for the pi-flux C4
bad = [round(x, 2) for x in np.arange(0.05, 1.0, 0.05) if (lambda E, O: E[1] < O[0] - 1e-9)(sorted(np.linalg.eigvalsh(ham(4, edges, x, [0, 0, 0, np.pi]))), sorted(np.linalg.eigvalsh(ham(4, edges, x, [0, 0, 0, np.pi])[np.ix_(odd, odd)])))]
print("s in 0.05..0.95 (step 0.05) where the pi-flux C4 violates Conjecture 2:", bad)
check("the violation holds for every s >= 0.6 on the grid", all(x in bad for x in [0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]))
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
