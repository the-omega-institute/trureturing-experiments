#!/usr/bin/env python3
"""Exact checks of the CMY singular fixed-sign LASSO intersection constant.

Only rational numbers and SymPy real algebraic expressions are used. These
checks do not establish the universal theorem or constitute a Lean probe.
"""
from itertools import combinations, product
from random import Random
import sympy as s

Q = s.Rational
rng = Random(261015241)


def simp(x):
    return s.simplify(x)


def sign(x):
    x = simp(x)
    if x == 0:
        return 0
    if x.is_positive is True:
        return 1
    if x.is_negative is True:
        return -1
    raise AssertionError(f"undecided exact algebraic sign: {x}")


def le(a, b):
    return sign(b - a) >= 0


def sq(v):
    return simp(v.dot(v))


def negsq(x):
    return sum(a*a for a in x if sign(a) < 0)


def subsets(items):
    items = tuple(items)
    for k in range(len(items) + 1):
        yield from combinations(items, k)


def vertices(A, p):
    b = A*p
    d = A.cols
    out = []
    for I in combinations(range(d), A.rank()):
        C = A[:, list(I)]
        if C.det() == 0:
            continue
        u = C.inv()*b
        if all(le(0, a) for a in u):
            v = s.zeros(d, 1)
            for i, a in zip(I, u):
                v[i] = a
            if v not in out:
                out.append(v)
    assert out
    return out


def faces(V):
    d = V[0].rows
    zero_sets = set()
    for vv in subsets(range(len(V))):
        if vv:
            zero_sets.add(tuple(j for j in range(d)
                                if all(V[i][j] == 0 for i in vv)))
    out = []
    for J in sorted(zero_sets):
        W = [v for v in V if all(v[j] == 0 for j in J)]
        p = sum(W, s.zeros(d, 1))/len(W)
        assert tuple(j for j in range(d) if p[j] == 0) == J
        out.append((J, p))
    return out


def spectrum(A, V):
    B = s.Matrix.hstack(*A.nullspace())
    H = B.T*B
    P = B*H.inv()*B.T
    cache = {}
    candidates = []
    for J, p in faces(V):
        for I in subsets(J):
            if not I:
                continue
            if I not in cache:
                D = s.diag(*[int(j in I) for j in range(A.cols)])
                E = H.inv()*B.T*D*B
                cache[I] = E.eigenvects()
            for eta, mult, basis in cache[I]:
                if sign(eta) <= 0 or sign(1-eta) <= 0:
                    continue
                # The selected cells have simple nontrivial eigenvalues.
                assert mult == 1, ('repeated eigenvalue requires a cone test', A, I)
                for eps in (1, -1):
                    v = (eps*B*basis[0]).applyfunc(simp)
                    if any((j in I) != (sign(v[j]) < 0) for j in J):
                        continue
                    if not all(le(v.dot(w-p), 0) for w in V):
                        continue
                    Dv = s.Matrix([v[j] if j in I else 0 for j in range(A.cols)])
                    assert (P*Dv-eta*v).applyfunc(simp) == s.zeros(A.cols, 1)
                    candidates.append((simp(eta), J, p, v, P))
    assert candidates, ('empty spectral candidate set', A, V)
    best = candidates[0]
    for c in candidates[1:]:
        if le(c[0], best[0]):
            best = c
    return best


def d_aff_sq(A, b, x):
    # A has independent rows in every test.
    r = A*x-b
    return simp((r.T*(A*A.T).inv()*r)[0])


def d_poly_sq(A, b, V, x):
    """Independent exact nearest-point computation by affine-face projection."""
    vals = []
    d = A.cols
    for J, _ in faces(V):
        C = A.col_join(s.eye(d)[list(J), :]) if J else A
        rhs = b.col_join(s.zeros(len(J), 1))
        # Drop dependent rows, retaining exactly their corresponding RHS.
        independent = list(C.T.rref()[1])
        C, rhs = C[independent, :], rhs[independent, :]
        z = x-C.T*(C*C.T).inv()*(C*x-rhs)
        if all(le(0, a) for a in z):
            assert A*z == b
            vals.append(sq(x-z))
    assert vals
    best = vals[0]
    for a in vals[1:]:
        if le(a, best):
            best = a
    return best


def ratio_bounded(r, eta):
    # r <= 2/(1-sqrt(1-eta)), using only order-preserving squaring.
    assert 0 < eta and eta < 1
    return le(r, 2) or le(eta*r*r-4*r+4, 0)


def attainment(A, V, best):
    eta, J, p, raw, P = best
    v = raw/s.sqrt(sq(raw))
    q0 = s.Matrix([max(0, a) if j in J else a for j, a in enumerate(v)])
    c = s.sqrt(1-eta)
    assert simp(sq(q0)-(1-eta)) == 0
    q = q0/c
    assert (P*q-c*v).applyfunc(simp) == s.zeros(A.cols, 1)
    r = v+q
    # Both v and q are unit vectors, so every coordinate of r has modulus <= 2.
    positive = [a for a in p if a > 0]
    delta = min(positive)/4 if positive else Q(1, 4)
    x = p+delta*r
    assert all(le(r.dot(w-p), 0) for w in V)
    assert all(sign(x[j]) == sign(r[j]) for j in J)
    assert all(le(0, x[j]) for j in range(A.cols) if j not in J)
    D = sq(x-p)
    a = d_aff_sq(A, A*p, x)
    b = simp(negsq(x))
    assert simp(a-b) == 0
    assert simp(a-delta*delta*eta) == 0
    assert simp(D/a-2*(1+c)/eta) == 0
    return simp(2*(1+c)/eta)


def fidelity():
    # Source Proposition 3: unique affine solution, feasible in the box/sign set.
    # n=2, A=2I, active={0}, lambda_n=1, y=(3,0).
    A = 2*s.eye(2)
    y = s.Matrix([3, 0])
    G = A.T*A/2
    M = s.diag(G[0, 0], 1)
    b = s.Matrix([(A.T*y/2)[0]-1, (A.T*y/2)[1]])
    p = M.inv()*b
    assert p == s.Matrix([1, 0])
    count = 0
    for a, g in product(range(-3, 4), repeat=2):
        x = s.Matrix([a, g])
        da = sq(x-p)
        dc = max(0, -a)**2+max(0, abs(g)-1)**2
        assert max(da, dc) == da
        count += 1
    assert sq(p) == 1  # x=0 attains ratio 1, so the least constant is 1.
    print(f'FIDELITY source Proposition 3, Eq. (20): tau=1; {count} exact points')


def family():
    e = s.symbols('e', positive=True)
    A = s.Matrix([[1, 1, 1], [1+e, 0, e]])
    v = s.Matrix([e, 1, -1-e])
    p = s.Matrix([0, 1, 1])
    y = A*p+s.Matrix([2, 0])
    assert (A*v).applyfunc(simp) == s.zeros(2, 1)
    assert (A.T*(A*p-y)/2+s.ones(3, 1)).applyfunc(simp) == s.zeros(3, 1)
    R = A*A.T
    assert simp(R.det()-(2+2*e+2*e*e)) == 0
    assert simp(s.trace(R)-(4+2*e+2*e*e)) == 0
    assert simp(8-s.trace(R)-2*(1-e)*(2+e)) == 0
    print('FAMILY symbolic: kernel=(e,1,-1-e); KKT exact; det=2+2e+2e^2, trace=4+2e+2e^2')
    print('FAMILY 0<e<=1: det>=2 and trace<=8 imply sigma_min(A)^2>=1/4')
    for e0 in [Q(1), Q(1, 2), Q(1, 3), Q(1, 5), Q(1, 10), Q(1, 100)]:
        C = A.subs(e, e0)
        V = vertices(C, p)
        best = spectrum(C, V)
        eta = e0*e0/(e0*e0+1+(1+e0)**2)
        assert simp(best[0]-eta) == 0
        t2 = attainment(C, V, best)
        assert ratio_bounded(4/(e0*e0), eta)
        assert le(4/(e0*e0), t2)
        print(f'FAMILY e={e0}, eta={eta}, tau^2={t2}, exact attaining point verified')


def small_cells():
    count = points = 0
    configs = []
    for _ in range(12):
        row = rng.sample(range(-5, 6), 3)
        configs.append((s.Matrix([[1]*3, row]), s.Matrix([rng.randrange(1, 4) for _ in row])))
    # Two-dimensional nullspaces, with interior and singleton feasible cells.
    for row in [[0, 1, 3, 7], [-4, -1, 2, 5], [-5, 0, 1, 4]]:
        configs.append((s.Matrix([[1]*4, row]), s.ones(4, 1)))
    configs.append((s.Matrix([[1]*4, [0, 1, 3, 7]]), s.zeros(4, 1)))
    for k, (A, p) in enumerate(configs):
        V = vertices(A, p)
        best = spectrum(A, V)
        eta = best[0]
        t2 = attainment(A, V, best)
        for _ in range(36):
            x = s.Matrix([Q(rng.randrange(-12, 17), 3) for _ in range(A.cols)])
            dp = d_poly_sq(A, A*p, V, x)
            da = d_aff_sq(A, A*p, x)
            dc = negsq(x)
            # Half the cases add one inactive zero design column and g in [-1,1].
            if k % 2:
                g = Q(rng.randrange(-6, 7), 2)
                dp += g*g
                da += g*g
                dc += max(0, abs(g)-1)**2
            den = max(da, dc)
            assert den != 0 or dp == 0
            if den:
                assert ratio_bounded(simp(dp/den), eta), (A, p, x, eta)
            points += 1
        print(f'CELL {k:02d}: n=2, active={A.cols}, vertices={len(V)}, eta={eta}, tau^2={t2}')
        count += 1
    print(f'SMALL_CELLS {count} deterministic random/fixed cells, {points} exact rational points; all attaining witnesses verified')


if __name__ == '__main__':
    fidelity()
    family()
    small_cells()
    print('ALL_CHECKS_OK')
