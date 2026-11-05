#!/usr/bin/env python3
"""Exact finite-atom checks for the higher even-moment tail-collapse candidate.

All quantities are Fractions.  The finite LP is solved by enumerating its rational
basic feasible solutions on a deliberately finite support; the global certificate
then identifies the same value for the unrestricted moment problem.
"""
from fractions import Fraction as F
from itertools import combinations


def solve_linear(A, b):
    """Return the unique rational solution, or None, by Gauss-Jordan elimination."""
    n, m = len(A), len(A[0]) if A else 0
    M = [[F(x) for x in row] + [F(y)] for row, y in zip(A, b)]
    row = 0
    pivots = []
    for col in range(m):
        pivot = next((i for i in range(row, n) if M[i][col]), None)
        if pivot is None:
            continue
        M[row], M[pivot] = M[pivot], M[row]
        q = M[row][col]
        M[row] = [x / q for x in M[row]]
        for i in range(n):
            if i != row and M[i][col]:
                q = M[i][col]
                M[i] = [u - q * v for u, v in zip(M[i], M[row])]
        pivots.append(col)
        row += 1
    if row < n:
        if any(all(M[i][j] == 0 for j in range(m)) and M[i][-1] != 0
               for i in range(row, n)):
            return None
        return None
    if len(pivots) != m:
        return None
    return [M[i][-1] for i in range(m)]


def finite_lp(support, m, K, t, two_sided):
    """Exact max over laws on support with E1=1,E X=0,E X²=1,E X^(2m)<=K."""
    best = None
    best_w = None
    # A vertex has 3 positive coordinates if the moment inequality is slack,
    # or 4 if it is active.  Smaller supports occur in these enumerations too.
    for active in (False, True):
        r = 4 if active else 3
        equations = [lambda x: F(1), lambda x: x, lambda x: x*x]
        rhs = [F(1), F(0), F(1)]
        if active:
            equations.append(lambda x: x ** (2*m))
            rhs.append(K)
        for inds in combinations(range(len(support)), r):
            pts = [support[i] for i in inds]
            A = [[fn(x) for x in pts] for fn in equations]
            w = solve_linear(A, rhs)
            if w is None or any(x < 0 for x in w):
                continue
            moment = sum(wi * x ** (2*m) for wi, x in zip(w, pts))
            if moment > K:
                continue
            hit = (lambda x: abs(x) >= t) if two_sided else (lambda x: x >= t)
            value = sum(wi for wi, x in zip(w, pts) if hit(x))
            if best is None or value > best:
                best, best_w = value, list(zip(pts, w))
    if best is None:
        raise AssertionError("finite LP found no feasible basic solution")
    return best, best_w


def case(m, t, a):
    T = t * t
    p = (1 - a) / (T - a)
    s2 = a
    # s is chosen as a rational square root in all cases below.
    s = {F(1, 4): F(1, 2), F(9, 16): F(3, 4), F(4, 9): F(2, 3)}[a]
    wplus = (1 - p - p * t / s) / 2
    wminus = (1 - p + p * t / s) / 2
    K = p * T**m + (1 - p) * a**m
    assert wplus >= 0 and wminus >= 0
    law = [(-s, wminus), (s, wplus), (t, p)]
    assert sum(w for _, w in law) == 1
    assert sum(x*w for x, w in law) == 0
    assert sum(x*x*w for x, w in law) == 1
    assert sum(x**(2*m)*w for x, w in law) == K
    # R_a(y)=y^m-m*a^(m-1)y+(m-1)a^m and Q=R_a(x^2)/R_a(T).
    D = T**m - m*a**(m-1)*T + (m-1)*a**m
    assert D > 0
    assert sum((x**(2*m) - m*a**(m-1)*x*x + (m-1)*a**m) * w
               for x, w in law) / D == p
    grid = set(F(x) for x in range(-4, 5))
    grid |= {F(-3, 4), F(-2, 3), F(-1, 2), F(-1, 4), F(1, 4), F(1, 2), F(2, 3), F(3, 4)}
    grid |= {t, -t, s, -s}
    grid = sorted(grid)
    one, one_w = finite_lp(grid, m, K, F(t), False)
    two, two_w = finite_lp(grid, m, K, F(t), True)
    assert one == p, (m, t, a, "one-sided", one, p, one_w)
    assert two == p, (m, t, a, "two-sided", two, p, two_w)
    # The finite LP checks also reproduce the source's m=2 tail formula.
    if m == 2:
        source = (K - 1) / ((T - 1)**2 + K - 1)
        assert source == p
    print(f"m={m} t={t} a={a}: K={K}, p={p}, LP1={one}, LP2={two}, support={law}")


def main():
    # Interior and boundary points in the far-tail range, for m=2,3,4.
    for m in (2, 3, 4):
        case(m, 2, F(9, 16))
        case(m, 3, F(1, 4))
        case(m, 3, F(4, 9))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
