#!/usr/bin/env python3
"""Exact checks for the set-theoretic question in ADK Example 4.6.

No finite test proves the all-n statement. Symbolic curves check the proposed
boundary mechanism; finite rational configurations check cofactors and ranks.
Requires Python 3 and SymPy. No floating-point arithmetic is used.
"""
from itertools import combinations, product
from random import Random
import sympy as S


def ensure(condition, message):
    if not condition:
        raise AssertionError(message)


def zero(expr):
    return S.cancel(expr) == 0


def same_projective(u, v):
    return any(not zero(a) for a in u) and any(not zero(a) for a in v) and all(
        zero(u[i] * v[j] - u[j] * v[i]) for i in range(3) for j in range(3)
    )


def segre(pairs):
    return S.Matrix([[x[a] * y[b] for a in range(3) for b in range(3)]
                     for x, y in pairs])


def cofactors(N):
    ensure(N.shape == (8, 9), 'cofactor dimensions')
    return S.Matrix([(-1) ** k * N[:, [j for j in range(9) if j != k]].det()
                     for k in range(9)])


def unvec(c):
    # Column vectorization: index 3*a+b is F[b,a].
    return S.Matrix(3, 3, lambda b, a: c[3 * a + b])


def equations(pairs):
    M = segre(pairs)
    minors9 = [M[list(J), :].det() for J in combinations(range(len(pairs)), 9)]
    cubes8 = [unvec(cofactors(M[list(I), :])).det()
              for I in combinations(range(len(pairs)), 8)]
    return all(zero(x) for x in minors9 + cubes8), len(minors9), len(cubes8)


A = S.eye(3).row_join(S.zeros(3, 1))
B = S.eye(3).row_join(S.Matrix([0, 0, 1]))
F0 = S.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
e = [S.eye(3)[:, i] for i in range(3)]
t = S.Symbol('t')

# Fidelity: the source's determinant block and its vectorized epipolar form.
x = S.Matrix(S.symbols('x0:3'))
y = S.Matrix(S.symbols('y0:3'))
block = A.row_join(x).row_join(S.zeros(3, 1)).col_join(
    B.row_join(S.zeros(3, 1)).row_join(y))
epipolar = (y.T * F0 * x)[0]
ensure(zero(block.det() + epipolar), 'source block determinant, global sign')
ensure(zero((segre([(x, y)]) * S.Matrix([F0[b, a] for a in range(3)
                                                       for b in range(3)]))[0]
            - epipolar), 'source vectorization')
print('FIDELITY: source fundamental block determinant and vectorization verified symbolically')

rng = Random(260631965)
pairs = []
while len(pairs) < 10:
    q = S.Matrix([rng.randint(-5, 5) for _ in range(4)])
    xx, yy = A * q, B * q
    if xx != S.zeros(3, 1) and yy != S.zeros(3, 1):
        pairs.append((xx, yy))
ensure(segre(pairs[:8]).rank() == 8, 'fidelity rank-eight premise')
c = cofactors(segre(pairs[:8]))
ensure(c != S.zeros(9, 1), 'nonzero source cofactor')
ensure(segre(pairs[:8]) * c == S.zeros(8, 1), 'source cofactor kernel')
ensure(unvec(c).rank() == 2, 'source rank-two eight-point matrix')
for n in [8, 9, 10]:
    ok, m9, c8 = equations(pairs[:n])
    ensure(ok, f'source camera equations n={n}')
    print(f'FIDELITY: actual camera tuple n={n}, rank=8, {m9} maximal minors and {c8} cofactor determinants vanish')

# A rank-eight tuple that fails the cubic: eight independent rows in x dot y=0.
nonimage = []
for xx in [S.Matrix(v) for v in product([-1, 0, 1], repeat=3) if any(v)]:
    for yy in [S.Matrix(v) for v in product([-1, 0, 1], repeat=3) if any(v)]:
        if (yy.T * xx)[0] == 0 and segre(nonimage + [(xx, yy)]).rank() > len(nonimage):
            nonimage.append((xx, yy))
            if len(nonimage) == 8:
                break
    if len(nonimage) == 8:
        break
ensure(len(nonimage) == 8, 'negative-control size')
ensure(unvec(cofactors(segre(nonimage))).det() != 0, 'negative control rejected by cubic')
ensure(not equations(nonimage + [nonimage[0]])[0], 'n=9 negative control')
print('CONTROL: rank-eight kernel spanned by I3 is rejected by cofactor determinant at n=8,9')

# Source rank-deficient boundary: the kernel consists only of rank-one matrices.
boundary = [(e[a], e[b]) for a in range(3) for b in range(3) if (a, b) != (0, 0)]
ensure(segre(boundary).rank() == 8, 'rank-one-only boundary rank')
ensure(unvec(segre(boundary).nullspace()[0]).rank() == 1, 'rank-one-only kernel')
for n in [9, 10]:
    pp = boundary + boundary[:n - 8]
    ok, m9, c8 = equations(pp)
    ensure(ok, 'rank-one boundary equations')
    print(f'BOUNDARY: n={n}, rank=8, kernel rank-one only; {m9} minors and {c8} cofactor determinants vanish')
lowrank = [(e[1], e[2])] * 9
ensure(segre(lowrank).rank() == 1 and equations(lowrank)[0], 'low-rank cofactor-zero branch')
print('BOUNDARY: n=9 rank-one Segre matrix has identically zero eight-row cofactors')

# Every rational projective pair in this finite grid compatible with E00.
points = []
for v in product([-1, 0, 1], repeat=3):
    if any(v) and next(a for a in v if a) == 1:
        points.append(S.Matrix(v))
ensure(len(points) == 13, 'projective test-grid size')
Ft = S.diag(1, t ** 2, 0)
branch_counts = [0, 0, 0]
for xx, yy in product(points, repeat=2):
    if xx[0] * yy[0] != 0:
        continue
    xt, yt = xx.copy(), yy.copy()
    if xx[0] != 0:
        branch = 0
        yt[0] = -t ** 2 * xx[1] * yy[1] / xx[0]
    elif yy[0] != 0:
        branch = 1
        xt[0] = -t ** 2 * xx[1] * yy[1] / yy[0]
    else:
        branch = 2
        xt[0], yt[0] = t * xx[1], -t * yy[1]
    branch_counts[branch] += 1
    ensure(xt.subs(t, 0) == xx and yt.subs(t, 0) == yy, 'boundary specialization')
    ensure(zero((yt.T * Ft * xt)[0]), 'rank-one deformation identity')
    for value in [S.Integer(-1), S.Integer(1), S.Integer(2)]:
        ensure(Ft.subs(t, value).rank() == 2, 'deformation matrix rank')
        ensure(xt.subs(t, value) != S.zeros(3, 1) and
               yt.subs(t, value) != S.zeros(3, 1), 'deformation nonzero images')
ensure(branch_counts == [36, 36, 16], 'all compatible grid pairs covered')
print(f'DEFORMATION: all 88 E00-compatible grid pairs; branch counts={branch_counts}; symbolic identities and nonzero specializations at t=-1,1,2')

# Check formulas with independent symbolic coordinates too.
a, b, c1, d, h = S.symbols('a b c d h', nonzero=True)
for xt, yt in [(S.Matrix([a, b, c1]), S.Matrix([-t*t*b*d/a, d, h])),
               (S.Matrix([-t*t*b*d/a, b, c1]), S.Matrix([a, d, h])),
               (S.Matrix([t*b, b, c1]), S.Matrix([-t*d, d, h]))]:
    ensure(zero((yt.T * Ft * xt)[0]), 'symbolic generic branch')
print('DEFORMATION: all three branches also verified with independent exact symbolic parameters')

# Rank-two reconstruction including every epipole branch on the same grid.
counts = [0, 0, 0, 0]
for xx, yy in product(points, repeat=2):
    if (yy.T * F0 * xx)[0] != 0:
        continue
    xepi = xx[0] == xx[1] == 0
    yepi = yy[0] == yy[1] == 0
    if not xepi and not yepi:
        counts[0] += 1
        k = 0 if xx[0] != 0 else 1
        lam = yy[k] / xx[k]
        ensure(lam != 0, 'reconstruction denominator')
        q = S.Matrix([xx[0], xx[1], xx[2], yy[2] / lam - xx[2]])
        ensure(same_projective(A*q, xx) and same_projective(B*q, yy), 'ordinary reconstruction')
    elif xepi and not yepi:
        counts[1] += 1
        xt = S.Matrix([t*yy[0], t*yy[1], xx[2]])
        q = S.Matrix([t*yy[0], t*yy[1], xx[2], t*yy[2]-xx[2]])
        ensure(A*q == xt and B*q == t*yy and xt.subs(t, 0) == xx, 'first epipole curve')
    elif yepi and not xepi:
        counts[2] += 1
        yt = S.Matrix([t*xx[0], t*xx[1], yy[2]])
        q = S.Matrix([t*xx[0], t*xx[1], t*xx[2], yy[2]-t*xx[2]])
        ensure(A*q == t*xx and B*q == yt and yt.subs(t, 0) == yy, 'second epipole curve')
    else:
        counts[3] += 1
        q = S.Matrix([0, 0, 1, 0])
        ensure(same_projective(A*q, xx) and same_projective(B*q, yy), 'both epipoles')
print(f'RECONSTRUCTION: all {sum(counts)} F0-compatible grid pairs; ordinary/first-epipole/second-epipole/both={counts}')
print('SCOPE: exact symbolic identities and bounded rational tests; no all-n proof or Zariski-closure certificate is claimed')
print('ALL_CHECKS_OK')
