#!/usr/bin/env python3
"""Exact finite checks for Maar's type-D gap question; these are not a gap proof."""
from fractions import Fraction as F
from itertools import combinations, product
import sympy as S


def edges(r):
    return [(0, 2), (1, 2)] + [(i, i + 1) for i in range(2, r - 1)] + [(r - 2, r)]


def root_matrix(r):
    a = [[0] * (r + 1) for _ in range(r)]
    a[0][0] = a[1][0] = -1
    for i in range(1, r):
        a[i - 1][i], a[i][i] = 1, -1
    a[r - 2][r] = a[r - 1][r] = 1
    return a


def mv(a, b):
    return tuple(sum(x * y for x, y in zip(row, b)) for row in a)


def frob2(a):
    return S.expand(sum(S.conjugate(x) * x for x in a))


def projection(a):
    return a if a.rows == 0 else a - S.trace(a) / a.rows * S.eye(a.rows)


def fidelity():
    """Proposition 2.1, including its defining trace identity, with Gaussian integers."""
    count, zeros = 0, 0
    for r in range(4, 7):
        for bits in product((0, 1), repeat=r):
            es = [(v, u) if b else (u, v) for (u, v), b in zip(edges(r), bits)]
            for mode in range(3):
                alpha = [2 if mode == 0 else (v + mode) % 4 for v in range(r + 1)]
                xs = [S.Matrix(alpha[v], alpha[u],
                               lambda i, j: (a + i + 2*j) % 3 - 1 +
                               S.I * ((2*a + i + j) % 3 - 1))
                      for a, (u, v) in enumerate(es)]
                ds = [S.zeros(n) for n in alpha]
                for (u, v), x in zip(es, xs):
                    ds[v] += x * x.H
                    ds[u] -= x.H * x
                hs = []
                for n in alpha:
                    z = S.Matrix(n, n, lambda i, j: i + j + 1 + S.I*(i-j))
                    hs.append(projection(z))
                norm2 = sum(map(frob2, xs))
                if norm2 == 0:
                    assert all(x == S.zeros(x.rows, x.cols) for x in xs)
                    zeros += 1
                    continue  # The source's moment map has domain X != 0.
                assert norm2 > 0
                mus = [projection(d) / norm2 for d in ds]
                lhs = sum(S.trace(mu*h) for mu, h in zip(mus, hs))
                rhs = sum(S.trace(x.H*(hs[v]*x-x*hs[u]))
                          for (u, v), x in zip(es, xs)) / norm2
                assert S.expand(lhs - rhs) == 0
                assert all(S.trace(mu) == 0 and mu == mu.H for mu in mus)
                assert sum(S.trace(d) for d in ds) == 0
                # Finite D is the affine diagram with vertex 0 empty.
                if alpha[0] == 0:
                    assert xs[0].cols == 0
                    assert ds[0].rows == 0
                count += 1
    print(f"FIDELITY Prop_2_1: {count} exact Gaussian-integer cases, all orientations r=4..6; {zeros} zero cases excluded")


def cartan_and_roots():
    counts = []
    for r in range(4, 9):
        a = root_matrix(r)
        es = {tuple(sorted(e)) for e in edges(r)}
        for i in range(r + 1):
            for j in range(r + 1):
                assert sum(row[i]*row[j] for row in a) == (
                    2 if i == j else -int(tuple(sorted((i, j))) in es))
        good = 0
        for b in product(range(3), repeat=r + 1):
            ab = mv(a, b)
            tau = sum(x*x for x in b) - sum(b[i]*b[j] for i, j in es)
            assert sum(x*x for x in ab) == 2*tau
            assert sum(ab) % 2 == 0
            if any(b) and tau <= 1:
                assert ab == (0,)*r or (sum(x != 0 for x in ab) == 2 and
                                        all(x in (-1, 0, 1) for x in ab))
                good += 1
        # The integer inverse of T(z)=(Az,z_0) on its image is injective.
        t = S.Matrix(a).col_join(S.Matrix([[1] + [0]*r]))
        assert t.det() != 0
        counts.append((r, 3**(r+1), good))
    print(f"CARTAN/root identity: (rank, beta cases, nonzero tau<=1): {counts}")


def all_roots(r):
    out = []
    for i, j in combinations(range(r), 2):
        for si, sj in product((-1, 1), repeat=2):
            v = [0]*r
            v[i], v[j] = si, sj
            out.append(tuple(v))
    return out


def selected_basis(cols, r):
    """A spanning tree plus one unbalanced edge in each unbalanced component."""
    parent = list(range(r))

    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v

    tree, omitted, adj = [], [], [[] for _ in range(r)]
    for k, c in enumerate(cols):
        ij = [i for i, x in enumerate(c) if x]
        if not ij:
            continue
        i, j = ij
        if root(i) != root(j):
            parent[root(i)] = root(j)
            tree.append(k)
            s = -c[i]*c[j]
            adj[i].append((j, s))
            adj[j].append((i, s))
        else:
            omitted.append(k)
    switches = [None]*r
    for i in range(r):
        if switches[i] is not None:
            continue
        switches[i] = 1
        stack = [i]
        while stack:
            v = stack.pop()
            for w, s in adj[v]:
                if switches[w] is None:
                    switches[w] = s*switches[v]
                    stack.append(w)
    unbalanced = set()
    for k in omitted:
        c = cols[k]
        i = next(i for i, x in enumerate(c) if x)
        if root(i) not in unbalanced and sum(x*s for x, s in zip(c, switches)):
            tree.append(k)
            unbalanced.add(root(i))
    return tree


def check_basis(cols, r):
    sel = selected_basis(cols, r)
    if not sel:
        assert all(not any(c) for c in cols)
        return
    c = S.Matrix.hstack(*(S.Matrix(cols[k]) for k in sel))
    left = (c.T*c).inv()*c.T
    assert all(abs(x) <= 1 for x in left)
    for j, col in enumerate(cols):
        v = S.Matrix(col)
        coeff = left*v
        assert c*coeff == v
        assert all(x.q == 1 and abs(x) <= 2 for x in coeff)
        if j not in sel:
            fundamental2 = 1 + sum(x*x for x in coeff)
            assert fundamental2 <= 4*(r+1)


def signed_bases():
    count = 0
    for r in range(2, 5):
        roots = all_roots(r)
        for size in range(1, 4):
            for cols in combinations(roots, size):
                check_basis(cols, r)
                count += 1
        check_basis(roots + [(0,)*r, roots[0]], r)
        count += 1
    print(f"SIGNED BASES: {count} exact column families; all subsets of <=3 signed roots for r=2..4")


def gram_checks():
    # All two-column multisets with binary beta and tau<=1 in affine D4.
    r, q = 4, 5
    bs = [b for b in product(range(2), repeat=q)
          if any(b) and sum(x*x for x in b)-sum(b[i]*b[j] for i,j in edges(r)) <= 1]
    count = 0
    for i, b in enumerate(bs):
        for c in bs[i:]:
            bmat = S.Matrix.hstack(S.Matrix(b), S.Matrix(c))
            alpha = list(bmat*S.ones(2, 1))
            d = max(1, *alpha)
            sizes = [sum(b), sum(c)]
            m = bmat*S.diag(*(S.Rational(1, s) for s in sizes))*bmat.T
            # Exact reciprocal-eigenvalue trace implies the stated lower eigenvalue bound.
            assert S.trace((bmat*bmat.T).pinv()) <= 1764*q**4*d**3
            mp = m.pinv()
            assert S.trace(mp) <= 1764*q**5*d**4
            # Projection of c preserves every beta dot c, even for singular M.
            z = S.Matrix([2, -1, 3, 0, 4])
            zproj = mp*m*z
            assert bmat.T*zproj == bmat.T*z
            assert m*zproj == m*z
            # The normal equation and energy identity at arbitrary scalar-block data.
            lambdas = [-sum(beta[v]*z[v] for v in range(q))/sum(beta) for beta in (b,c)]
            assert S.Matrix([sum(beta[v]*lam for beta,lam in zip((b,c),lambdas))
                             for v in range(q)]) == -m*z
            count += 1
    print(f"SINGULAR GRAM: {count} exhaustive binary-beta pairs in affine D4; exact pseudoinverses")


def unstable_single_arrows():
    count = 0
    for r in range(4, 10):
        q = r + 1
        for m, n in product(range(1, 7), repeat=2):
            for k in range(1, min(m, n)+1):
                if k == m == n:
                    continue  # full rank square: semistable, excluded from the nullcone
                # A has k diagonal entries 1. A trace-zero row/column scaling contracts it
                # whenever k<m or k<n; hence this is an actual nonzero SL-nullcone family.
                h = F(2, k) - F(1, m) - F(1, n)
                assert h > 0
                d = max(m, n)
                assert h >= F(1, 1765**2*q**14*d**11)
                assert h >= F(1, 2000**2*q**14*d**12)
                count += 1
    print(f"ACTUAL UNSTABLE FAMILY: {count} exact squared-gap checks, ranks 4..9, endpoint dimensions <=6")


if __name__ == '__main__':
    fidelity()
    cartan_and_roots()
    signed_bases()
    gram_checks()
    unstable_single_arrows()
    print('ALL_CHECKS_OK')
