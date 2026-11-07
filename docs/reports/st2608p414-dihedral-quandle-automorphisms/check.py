#!/usr/bin/env python3
"""Independent exact checks of the proposed answer to Sangare--Ta Problem 4.14.

All coefficients are integers. Finite tests check the formulas, not exhaustion
of the infinite automorphism group and not a universal proof.
"""
from itertools import permutations, product
from math import gcd
from random import Random


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def mm(a, b):
    bt = list(zip(*b))
    return [[sum(x * y for x, y in zip(row, col)) for col in bt] for row in a]


def mv(a, v):
    return [sum(x * y for x, y in zip(row, v)) for row in a]


def det(a):
    """Bareiss elimination, with every division checked to be exact."""
    a = [row[:] for row in a]
    n, sign, old = len(a), 1, 1
    for k in range(n - 1):
        p = next((i for i in range(k, n) if a[i][k]), None)
        if p is None:
            return 0
        if p != k:
            a[p], a[k], sign = a[k], a[p], -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % old == 0
                a[i][j] = numerator // old
            a[i][k] = 0
        old = pivot
    return sign * a[-1][-1] if n else 1


def signed(m, k):
    n, N = 2 * m, 4 * m
    k %= N
    v = [0] * n
    v[k % n] = 1 if k < n else -1
    return v


def columns(vs):
    return list(map(list, zip(*vs)))


def affine(N, g, h):
    u, v = g
    a, b = h
    return u * a % N, (u * b + v) % N


def ainv(N, g):
    u, v = g
    a = pow(u, -1, N)
    return a, -a * v % N


def P(m, g):
    u, v = g
    return columns([signed(m, u * j + v) for j in range(2 * m)])


def A(m, g):
    N = 4 * m
    u, v = g
    return [[int(i == (u * j + v) % N) for j in range(N)] for i in range(N)]


def in_order(m, M):
    n = 2 * m
    T = columns([signed(m, j + 2) for j in range(n)])
    J = columns([signed(m, -j) for j in range(n)])
    return mm(M, T) == mm(T, M) and mm(M, J) == mm(J, M)


def in_K(m, M):
    n = 2 * m
    return (in_order(m, M) and abs(det(M)) == 1
            and all((M[i][j] - int(i == j)) % 2 == 0
                    for i in range(n) for j in range(n)))


def lift(m, M):
    n, N = 2 * m, 4 * m
    out = eye(N)
    for j in range(N):
        s, sign = j % n, 1 if j < n else -1
        for r in range(n):
            x = sign * (M[r][s] - int(r == s))
            assert x % 2 == 0
            out[r][j] += x // 2
            out[r + n][j] -= x // 2
    return out


def algebra_product(x, y):
    N = len(x)
    out = [0] * N
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                if b:
                    out[(2 * j - i) % N] += a * b
    return out


def preserves_product(M):
    N = len(M)
    cols = list(map(list, zip(*M)))
    return all(algebra_product(cols[i], cols[j]) == cols[(2 * j - i) % N]
               for i in range(N) for j in range(N))


def shear(m, k, upper):
    n = 2 * m
    out = eye(n)
    for j in range(n):
        if j % 2 == int(upper):
            a, b = signed(m, j - 1), signed(m, j + 1)
            for r in range(n):
                out[r][j] += 2 * k * (a[r] + b[r])
    return out


def fidelity_R4():
    m, N = 1, 4
    perms = []
    for p in permutations(range(N)):
        M = [[int(i == p[j]) for j in range(N)] for i in range(N)]
        if preserves_product(M):
            perms.append(M)
    a = [[int(i == [0, 3, 2, 1][j]) for j in range(N)] for i in range(N)]
    b = [[int(i == [2, 1, 0, 3][j]) for j in range(N)] for i in range(N)]
    tau = [[int(i == [1, 0, 3, 2][j]) for j in range(N)] for i in range(N)]
    assert len(perms) == 8
    assert mm(a, a) == mm(b, b) == mm(tau, tau) == eye(N)
    assert mm(a, b) == mm(b, a) and mm(mm(tau, a), tau) == b
    K = [[[x, 0], [0, y]] for x, y in product([-1, 1], repeat=2)]
    F = [(u, v) for u in [1, 3] for v in range(N)]
    images = {tuple(map(tuple, mm(lift(m, M), A(m, g)))) for M in K for g in F}
    assert images == {tuple(map(tuple, M)) for M in perms}
    assert len(K) * len(F) // 4 == 8
    print('FIDELITY R4: 8 permutation automorphisms, BPS generators and D4 relations; model images=8')


def check_m(m):
    n, N = 2 * m, 4 * m
    F = [(u, v) for u in range(N) if gcd(u, N) == 1 for v in range(N)]
    L = [(u, v) for u in [1, 1 + n] for v in [0, n]]
    ident = eye(n)
    rng = Random(414000 + m)
    pairs = [((1, 0), (1, 1)), (F[-1], F[len(F) // 2])]
    pairs += [(rng.choice(F), rng.choice(F)) for _ in range(12)]
    for g, h in pairs:
        assert mm(P(m, g), P(m, h)) == P(m, affine(N, g, h))
        assert mm(A(m, g), A(m, h)) == A(m, affine(N, g, h))
        assert mm(P(m, g), P(m, ainv(N, g))) == ident
    # Exhaust every affine parameter, including the exact four-element quotient kernel.
    reduction_kernel = [g for g in F if g[0] % n == 1 % n and g[1] % n == 0]
    assert reduction_kernel == L
    for g in F:
        assert in_order(m, mm(mm(P(m, g), shear(m, 1, True)), P(m, ainv(N, g))))
    for l in L:
        pinv = P(m, ainv(N, l))
        assert in_K(m, pinv)
        assert mm(lift(m, pinv), A(m, l)) == eye(N)
    samples = []
    for k in [-2, -1, 1, 3]:
        for upper in [False, True]:
            M, inv = shear(m, k, upper), shear(m, -k, upper)
            assert in_K(m, M) and mm(M, inv) == ident
            assert mm(lift(m, M), lift(m, inv)) == eye(N)
            assert preserves_product(lift(m, M))
            samples.append(M)
    for _ in range(6):
        M, H = rng.choice(samples), rng.choice(samples)
        g, h = rng.choice(F), rng.choice(F)
        conj = mm(mm(P(m, g), H), P(m, ainv(N, g)))
        combined = mm(M, conj)
        assert in_K(m, combined)
        left = mm(mm(lift(m, M), A(m, g)), mm(lift(m, H), A(m, h)))
        right = mm(lift(m, combined), A(m, affine(N, g, h)))
        assert left == right and preserves_product(left)
    # Check a spanning family of the integral commutant equations in the R^2 description.
    # A signed cyclic block acts by multiplication in Z[t]/(t^m+1).
    def shift_poly(p, power):
        out = [0] * m
        for i, x in enumerate(p):
            q, r = divmod(i + power, m)
            out[r] += (-1 if q % 2 else 1) * x
        return out
    def star(p):
        out = [0] * m
        for i, x in enumerate(p):
            q, r = divmod(-i, m)
            out[r] += (-1 if q % 2 else 1) * x
        return out
    def block(polys):
        out = [[0] * n for _ in range(n)]
        for r in range(2):
            for s in range(2):
                for j in range(m):
                    col = shift_poly(polys[2*r+s], j)
                    for i, x in enumerate(col):out[2*i+r][2*j+s] = x
        return out
    for j in range(m):
        z, e = [0] * m, [int(i == j) for i in range(m)]
        for slot in range(4):
            p = [x[:] for x in [z,z,z,z]]
            st = star(e)
            if slot == 1: st = shift_poly(st, 1)
            if slot == 2: st = shift_poly(st, -1)
            p[slot] = [x+y for x,y in zip(e,st)]
            assert in_order(m, block(p))
    for _ in range(40):
        ps = [[rng.randrange(-1, 2) for _ in range(m)] for _ in range(4)]
        equations = (ps[0] == star(ps[0]) and ps[3] == star(ps[3])
                     and ps[1] == shift_poly(star(ps[1]), 1)
                     and ps[2] == shift_poly(star(ps[2]), -1))
        assert in_order(m, block(ps)) == equations
    # Directly form S_f, square it, and compute the diagonal trace, over integers.
    d = [rng.randrange(-3, 4) for _ in range(n)]
    S = [[0] * N for _ in range(N)]
    for j, x in enumerate(d):
        for i in range(N):S[(2*j-i) % N][i] += x
    trace_square = sum(S[i][j]*S[j][i] for i in range(N) for j in range(N))
    assert trace_square == N * sum(x*x for x in d)
    if m >= 2:
        assert shear(m, 1, True) != ident
        assert mm(shear(m, 2, True), shear(m, 3, True)) == shear(m, 5, True)
    print(f'm={m}: affine={len(F)}, kernel={len(L)}, lifts=8, mixed_products=6, '
          f'basis_pairs_per_lift={N*N}, order_blocks={4*m}, arbitrary_blocks=40, trace={trace_square}; OK')


def main():
    fidelity_R4()
    for m in range(1, 9):
        check_m(m)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
