#!/usr/bin/env python3
"""Exact finite checks for Wallman--Flammia, arXiv:1404.6025v4, Section VII.

Python 3 + SymPy. All evidence uses rationals and exact algebraic numbers.
The finite samples support the proposed argument; they do not prove the conjecture.
"""
from itertools import product
from functools import cache
import sympy as S
from sympy.polys.matrices import DomainMatrix
from sympy.polys.densearith import dup_mul

R, I, z = S.Rational, S.I, S.Symbol("z")
omega = -R(1, 2) + I * S.sqrt(3) / 2
K = S.QQ.algebraic_field(omega)
ZERO, ONE = K.zero, K.one


@cache
def field_conjugate(v):
    return K.from_sympy(S.conjugate(K.to_sympy(v)))


def simp(x):
    return S.simplify(S.expand_complex(x))


def equal(A, B):
    assert A.shape == B.shape
    assert all(simp(x) == 0 for x in A - B)


def vec(A):
    return S.Matrix([A[r, c] for c in range(A.cols) for r in range(A.rows)])


def unit(n, i, j):
    A = S.zeros(n)
    A[i, j] = 1
    return A


def swap(d):
    return S.Matrix(d*d, d*d,
                    lambda r, c: int(r // d == c % d and r % d == c // d))


def haar_two_twirl(A, d):
    F = swap(d)
    a, b = S.trace(A), S.trace(F*A)
    return ((a-b/d)*S.eye(d*d)+(b-a/d)*F)/(d*d-1)


def kraus_map(ks, A):
    return sum((k*A*k.H for k in ks), S.zeros(A.rows)).applyfunc(simp)


def transfer(ks, bs):
    d = bs[0].rows
    equal(sum((k.H*k for k in ks), S.zeros(d)), S.eye(d))
    return S.Matrix(len(bs), len(bs),
                    lambda a, b: simp(S.trace(bs[a].H*kraus_map(ks, bs[b]))/d))


def monomial_action(U, bs):
    L = S.Matrix(len(bs), len(bs),
                 lambda a, b: simp(S.trace(bs[a].H*U*bs[b]*U.H)/U.rows))
    out = []
    for c in range(L.cols):
        nz = [(r, K.from_sympy(L[r, c])) for r in range(L.rows) if L[r, c] != 0]
        assert len(nz) == 1
        out.append(nz[0])
    return tuple(out)


def compose(a, b):
    return tuple((a[j][0], v*a[j][1]) for j, v in b)


def generated(gens):
    ident = tuple((j, ONE) for j in range(len(gens[0])))
    seen, todo = {ident}, [ident]
    while todo:
        a = todo.pop()
        for b in gens:
            c = compose(a, b)
            if c not in seen:
                seen.add(c)
                todo.append(c)
    return sorted(seen, key=repr)


def twirl_channel_square(L, group):
    """Direct uniform average of (O L O^-1) tensor (O L O^-1)."""
    n = L.rows
    nz = [(r, c, K.from_sympy(L[r, c]))
          for r in range(n) for c in range(n) if L[r, c] != 0]
    Q = [[ZERO for _ in range(n*n)] for _ in range(n*n)]
    gram = [[ZERO for _ in range(n-1)] for _ in range(n-1)]
    for g in group:
        ag = [(g[r][0], g[c][0], v*g[r][1]/g[c][1]) for r, c, v in nz]
        for r, c, v in ag:
            for a, b, w in ag:
                Q[n*r+a][n*c+b] += v*w
                if r == a and r and c and b:
                    gram[c-1][b-1] += field_conjugate(v)*w
    size = K.convert(len(group))
    return [[v/size for v in row] for row in Q], [[v/size for v in row] for row in gram]


def schur_stable(coeffs):
    """Exact real Schur recursion: all polynomial roots are strictly in |z|<1."""
    a = [S.Rational(x) for x in coeffs]
    while len(a) > 1:
        leading, constant = a[0], a[-1]
        assert abs(constant) < abs(leading), a
        a = [leading*a[j]-constant*a[-1-j] for j in range(len(a)-1)]
        a = [x/a[0] for x in a]


def test_channel(name, ks, bs, labels, group, unitary=False):
    d, n = bs[0].rows, len(bs)
    L = transfer(ks, bs)
    equal(L[:1, :], S.Matrix([[1]+[0]*(n-1)]))
    A = L[1:, 1:]
    f = simp(S.trace(A)/(n-1))
    u = simp(S.trace(A.H*A)/(n-1))
    assert u.is_Rational and 0 <= u <= 1
    J = S.Matrix(d*d, d*d,
                 lambda r, c: kraus_map(ks, unit(d, r//d, c//d))[r % d, c % d]/d)
    tau = kraus_map(ks, S.eye(d))/d
    assert simp((n-1)*u-d*d*S.trace(J*J)+d*S.trace(tau*tau)) == 0
    Q, gram = twirl_channel_square(L, group)
    assert gram == [[K.from_sympy(u) if i == j else ZERO
                     for j in range(n-1)] for i in range(n-1)]
    sectors = {}
    for a, x in enumerate(labels):
        for b, y in enumerate(labels):
            charge = tuple((x[j]+y[j]) % d for j in range(2))
            sectors.setdefault(charge, []).append(n*a+b)
    charge_of = {i: q for q, inds in sectors.items() for i in inds}
    assert all(Q[i][j] == ZERO for i in range(n*n) for j in range(n*n)
               if charge_of[i] != charge_of[j])
    coeffs = [ONE]
    for inds in sectors.values():
        block = DomainMatrix.from_list([[Q[i][j] for j in inds] for i in inds], K)
        coeffs = dup_mul(coeffs, block.charpoly(), K)
    coeffs = [simp(K.to_sympy(c)) for c in coeffs]
    assert all(c.is_Rational for c in coeffs)
    p = S.Poly.from_list(coeffs, z)
    residual, rem = S.div(p, S.Poly(z-1, z))
    assert rem.is_zero
    if unitary:
        assert u == 1 and residual.eval(1) == 0
    else:
        assert u < 1 and f*f <= u and residual.eval(1) != 0
        # Factoring over QQ preserves every root and its algebraic multiplicity.
        # Apply the exact criterion to each factor to control coefficient growth.
        _, factors = residual.factor_list()
        assert sum(poly.degree()*mult for poly, mult in factors) == residual.degree()
        for poly, mult in factors:
            assert mult > 0
            schur_stable(poly.all_coeffs())
    # Independently check the proposed diagonal block factorization.
    mm = [n*a+b for a in range(1, n) for b in range(1, n)]
    mcoeffs = [ONE]
    for inds0 in sectors.values():
        inds = [i for i in inds0 if i in mm]
        if inds:
            mcoeffs = dup_mul(mcoeffs, DomainMatrix.from_list(
                [[Q[i][j] for j in inds] for i in inds], K).charpoly(), K)
    pm = S.Poly.from_list([simp(K.to_sympy(c)) for c in mcoeffs], z)
    assert p == S.Poly((z-1)*(z-f)**(2*(n-1))*pm.as_expr(), z)
    if d == 2:
        # Source Lemma 6, Eq. (31): average gate fidelity and normalized Choi overlap.
        states = [(S.eye(2)+sign*B)/2 for B in bs[1:] for sign in [-1, 1]]
        favg = sum(S.trace(rho*kraus_map(ks, rho)) for rho in states)/6
        phi = vec(S.eye(2))/S.sqrt(2)
        ent = (phi.H*J*phi)[0]
        assert simp(3*favg-2*ent-1) == 0
        # Entire 3x3 matrix-space vectorization, independently from Kronecker entries.
        for r, c in product(range(3), repeat=2):
            E = unit(3, r, c)
            T = S.zeros(3)
            for g in group:
                O = S.zeros(3)
                for j in range(1, 4):
                    O[g[j][0]-1, j-1] = K.to_sympy(g[j][1])
                Ag = O*A*O.T
                T += Ag*E*Ag.T/len(group)
            mv = S.Matrix([[K.to_sympy(Q[4*a+b][4*(r+1)+c+1])]
                           for a in range(1, 4) for b in range(1, 4)])
            equal(mv, S.Matrix([T[a, b] for a in range(3) for b in range(3)]))
    print(f"PASS {name}: d={d}, |adjoint group|={len(group)}, f={f}, u={u}, "
          f"{'unitary boundary: repeated 1' if unitary else 'simple 1; remaining roots Schur-stable'}")


def main():
    paulis = [S.eye(2), S.Matrix([[0, 1], [1, 0]]),
              S.Matrix([[0, -I], [I, 0]]), S.diag(1, -1)]
    h = S.Matrix([[1, 1], [1, -1]])/S.sqrt(2)
    phase = S.diag(1, I)
    g2 = generated([monomial_action(U, paulis) for U in [h, phase]])
    assert len(g2) == 24
    # Full 2-design certification, on all 16 physical two-copy matrix units.
    # Clifford lifts of h and phase generate a finite group; global phases cancel.
    lifts = {tuple((j, ONE) for j in range(4)): S.eye(2)}
    todo = list(lifts)
    gens2 = [(monomial_action(U, paulis), U) for U in [h, phase]]
    while todo:
        a = todo.pop()
        for b, U in gens2:
            c = compose(b, a)
            if c not in lifts:
                lifts[c] = (U*lifts[a]).applyfunc(simp)
                todo.append(c)
    mom = S.zeros(16)
    for U in lifts.values():
        UU = S.kronecker_product(U, U).applyfunc(simp)
        mom += S.kronecker_product(UU.conjugate(), UU)/24
    mom = mom.applyfunc(simp)
    for a, b in product(range(4), repeat=2):
        E = unit(4, a, b)
        equal(mom*vec(E), vec(haar_two_twirl(E, 2)))
    print("PASS qubit exact 2-design: all 16 two-copy matrix units")
    labels2 = [(0, 0), (1, 0), (1, 1), (0, 1)]
    for p in [R(0), R(1, 4), R(1, 2), R(3, 4), R(1)]:
        test_channel(f"qubit phase-flip p={p}",
                     [S.sqrt(1-p)*paulis[0], S.sqrt(p)*paulis[3]],
                     paulis, labels2, g2, unitary=p in [0, 1])
    for s, t in [(0, 1), (R(3, 5), R(4, 5)), (R(4, 5), R(3, 5)), (1, 0)]:
        s, t = R(s), R(t)
        test_channel(f"qubit damping s={s}", [S.diag(1, s), t*unit(2, 0, 1)],
                     paulis, labels2, g2, unitary=s == 1)

    labels3 = list(product(range(3), repeat=2))
    X = S.Matrix(3, 3, lambda r, c: int(r == (c+1) % 3))
    Z = S.diag(1, omega, omega**2).applyfunc(simp)
    ws = [(omega**(2*a*b)*X**a*Z**b).applyfunc(simp) for a, b in labels3]
    F = S.Matrix(3, 3, lambda r, c: omega**(r*c)/S.sqrt(3)).applyfunc(simp)
    P = S.diag(1, omega**2, omega**2).applyfunc(simp)
    for U in [X, Z, F, P]:
        equal(U.H*U, S.eye(3))
        assert simp(S.det(U)**12) == 1
    g3 = generated([monomial_action(U, ws) for U in [X, Z, F, P]])
    assert len(g3) == 216
    # Translation phases kill nonzero total Weyl charge; SL(2,3) is transitive
    # on the eight nonzero labels. Verify this directly for every basis pair.
    for a, b in product(range(9), repeat=2):
        coeff = {}
        for g in g3:
            key = (g[a][0], g[b][0])
            coeff[key] = coeff.get(key, ZERO)+g[a][1]*g[b][1]/K.convert(216)
        avg = S.zeros(9)
        for (r, c), v in coeff.items():
            if v:
                avg += K.to_sympy(v)*S.kronecker_product(ws[r], ws[c])
        equal(avg, haar_two_twirl(S.kronecker_product(ws[a], ws[b]), 3))
    print("PASS qutrit exact 2-design: all 81 Weyl tensor basis elements; 216 adjoint actions")
    for s, t in [(0, 1), (R(3, 5), R(4, 5)), (R(4, 5), R(3, 5))]:
        s, t = R(s), R(t)
        ks = [S.diag(1, s, 1), t*unit(3, 0, 1)]
        # The unital projection is not even positive: its (0,0) output on |2><2| is negative.
        shift = kraus_map(ks, S.eye(3))/3-S.eye(3)/3
        assert simp((kraus_map(ks, unit(3, 2, 2))-shift)[0, 0]) == -t*t/3
        assert -t*t/3 < 0
        test_channel(f"qutrit selective damping s={s}", ks, ws, labels3, g3)
        test_channel(f"qutrit Fourier-after-damping s={s}",
                     [(F*k).applyfunc(simp) for k in ks], ws, labels3, g3)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
