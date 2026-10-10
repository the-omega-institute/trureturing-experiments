#!/usr/bin/env python3
"""Cruz Ortega, Menendez-Conde, arXiv:2306.08187v1, abstract and Section 5: "for every graph the maximum
[of lambda_2(G; w) subject to R(G; w) = R(G)] is attained for weights that are invariant under automorphisms".
R is the total effective resistance (Kirchhoff index, sum over all unordered vertex pairs).

Counterexample family K_{2,m}, m >= 3 (hubs u, v; leaves x_1..x_m). With conductance 2 on u-x_i and 1 on v-x_i
(Laplacian B_m), lambda_2(B_m) = 3 and R(B_m) = (2m^3+2m^2+5m+9)/(6m); the unweighted graph has lambda_2 = 2 and
R = D_m/(2m), D_m = m^3+m^2+2m+4. Scaling B_m by c_m = (2m^3+2m^2+5m+9)/(3 D_m) meets the resistance budget with
lambda_2 = 2 + (m+1)/D_m > 2, while the automorphism group is transitive on edges, so every invariant feasible
assignment is the unweighted one (lambda_2 = 2).

Checks (exact, sympy): for m = 3..8 the Laplacian spectra, the pseudoinverse Kirchhoff index, the harmonic-current
resistances of the source's definition for every pair class, the scaled assignment, and edge-transitivity of the
automorphism group (m = 3, 4 by enumeration)."""
import sys, itertools
import sympy as sp

bad = 0
def report(label, ok):
    global bad
    print(label, ok)
    bad += 0 if ok else 1

def laplacian(m, wu, wv):
    n = m + 2; L = sp.zeros(n, n)
    for i in range(m):
        for (h, w) in ((0, wu), (1, wv)):
            x = 2 + i
            L[h, h] += w; L[x, x] += w; L[h, x] -= w; L[x, h] -= w
    return L

def eff_res_harmonic(L, i, j):
    """Source definition: voltage 1 at i, 0 at j, harmonic elsewhere; r = 1 / (current out of i)."""
    n = L.shape[0]; interior = [k for k in range(n) if k not in (i, j)]
    vsym = sp.symbols(f'h0:{n}')
    sub = {vsym[i]: 1, vsym[j]: 0}
    eqs = [sum(L[k, l] * vsym[l] for l in range(n)) for k in interior]
    sol = sp.solve([e.subs(sub) for e in eqs], [vsym[k] for k in interior], dict=True)[0]
    sub.update(sol)
    current = sum(L[i, l] * sp.sympify(vsym[l]).subs(sub) for l in range(n))
    return sp.nsimplify(1 / current)

def kirchhoff_pinv(L):
    n = L.shape[0]; J = sp.ones(n, n) / n
    Lp = (L + J).inv() - J
    return sp.nsimplify(n * Lp.trace())

for m in range(3, 9):
    n = m + 2
    Dm = m**3 + m**2 + 2 * m + 4
    B = laplacian(m, 2, 1); L0 = laplacian(m, 1, 1)
    ev_B = sorted(sp.Matrix(B).eigenvals(multiple=True), key=lambda e: float(e))
    ev_0 = sorted(sp.Matrix(L0).eigenvals(multiple=True), key=lambda e: float(e))
    ok = sp.simplify(ev_B[1] - 3) == 0 and sp.simplify(ev_0[1] - 2) == 0
    RB = kirchhoff_pinv(B); R0 = kirchhoff_pinv(L0)
    ok &= RB == sp.Rational(2 * m**3 + 2 * m**2 + 5 * m + 9, 6 * m) and R0 == sp.Rational(Dm, 2 * m)
    # source harmonic definition, one representative of each pair class
    reps = {(0, 1): 1, (0, 2): m, (1, 2): m, (2, 3): m * (m - 1) // 2}
    RB_h = sum(cnt * eff_res_harmonic(B, a, b) for (a, b), cnt in reps.items())
    R0_h = sum(cnt * eff_res_harmonic(L0, a, b) for (a, b), cnt in reps.items())
    ok &= RB_h == RB and R0_h == R0
    c = sp.Rational(2 * m**3 + 2 * m**2 + 5 * m + 9, 3 * Dm)
    W = laplacian(m, 2 * c, c)
    RW = kirchhoff_pinv(W)
    lamW = sorted(sp.Matrix(W).eigenvals(multiple=True), key=lambda e: float(e))[1]
    ok &= RW == R0 and sp.simplify(lamW - (2 + sp.Rational(m + 1, Dm))) == 0 and lamW > 2
    print(f'm={m}: lambda2(B)={ev_B[1]}, R(B)={RB}, R(L0)={R0}, scaled lambda2={sp.nsimplify(lamW)} > 2, R(w)=R(L0): {RW == R0}')
    report(f'm={m}: all identities', ok)

# automorphism edge-transitivity
for m in (3, 4):
    n = m + 2
    edges = {frozenset((h, 2 + i)) for h in (0, 1) for i in range(m)}
    orbit_ok = True
    autos = [p for p in itertools.permutations(range(n)) if {frozenset((p[a], p[b])) for a, b in map(tuple, edges)} == edges]
    e0 = frozenset((0, 2))
    orbit = {frozenset((p[a], p[b])) for p in autos for (a, b) in [tuple(e0)]}
    report(f'm={m}: {len(autos)} automorphisms act transitively on the {len(edges)} edges:', orbit == edges)
print('bad=', bad)
sys.exit(1 if bad else 0)
