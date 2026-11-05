#!/usr/bin/env python3
"""Exact finite diagnostics for the positive exponential Krawtchouk radius.

The reference coefficient formula is evaluated independently of the matrix recurrence.
Perturbations are deterministic samples, not an exhaustive test of a real coefficient ball.
Only Python integers and fractions.Fraction participate in the checks.
"""
from fractions import Fraction as Q
from math import comb, factorial


def eye(size):
    return [[Q(i == j) for j in range(size)] for i in range(size)]


def sub(A, B):
    return [[x - y for x, y in zip(a, b)] for a, b in zip(A, B)]


def scale(q, A):
    return [[q * x for x in row] for row in A]


def norm(A):
    return max(sum(abs(x) for x in row) for row in A)


def jacobi(v):
    size = len(v)
    return [[v[i][0] if j == i + 1 else v[i][1] if j == i
             else v[i][2] if i == j + 1 else Q(0)
             for j in range(size)] for i in range(size)]


def left_mul(v, A):
    size = len(v)
    return [[v[i][1] * A[i][j]
             + (v[i][0] * A[i + 1][j] if i + 1 < size else 0)
             + (v[i][2] * A[i - 1][j] if i > 0 else 0)
             for j in range(size)] for i in range(size)]


def matrices(v):
    N = len(v) - 1
    B = [eye(N + 1), jacobi(v)]
    for m in range(1, N):
        a, b, c = v[m]
        JB = left_mul(v, B[m])
        B.append([[((JB[i][j] - b * B[m][i][j]) - c * B[m - 1][i][j]) / a
                   for j in range(N + 1)] for i in range(N + 1)])
    return B


def reference(N, r):
    return [(Q(N - i, N), (1 - r) * Q(i, N), r * Q(i, N))
            for i in range(N + 1)]


def admissible(v):
    N = len(v) - 1
    assert N >= 1 and v[0] == (1, 0, 0) and v[N][0] == 0
    assert all(sum(row) == 1 and min(row) >= 0 for row in v)
    assert all(v[i][0] > 0 for i in range(N))
    assert all(v[i][2] > 0 for i in range(1, N + 1))


def triangle(m, n, k, N):
    return abs(m - n) <= k <= min(m + n, N)


def tensor_coefficient(N, r, m, n, k):
    """Average Y_S Y_T using Y_i^2 = r + (1-r)Y_i (P1, cor:kraw)."""
    result = Q(0)
    for h in range(min(m, n) + 1):
        ell = k - m - n + 2 * h
        if 0 <= n - h <= N - m and 0 <= ell <= h:
            result += (Q(comb(m, h) * comb(N - m, n - h) * comb(h, ell), comb(N, n))
                       * (1 - r) ** ell * r ** (h - ell))
    return result


def samples(v0, eta):
    """Each one-row hexagon vertex, both terminal vertices, and six coupled arrays."""
    N = len(v0) - 1
    directions = [(1, -1, 0), (1, 0, -1), (0, 1, -1),
                  (-1, 1, 0), (-1, 0, 1), (0, -1, 1)]
    yield list(v0)
    for i in range(1, N):
        for direction in directions:
            v = list(v0)
            v[i] = tuple(x + eta * d for x, d in zip(v0[i], direction))
            yield v
    for sign in [-1, 1]:
        v = list(v0)
        v[N] = (Q(0), v0[N][1] + sign * eta, v0[N][2] - sign * eta)
        yield v
    for shift in range(6):
        v = list(v0)
        for i in range(1, N):
            direction = directions[(i + shift) % 6]
            v[i] = tuple(x + eta * d for x, d in zip(v0[i], direction))
        sign = (-1) ** shift
        v[N] = (Q(0), v0[N][1] + sign * eta, v0[N][2] - sign * eta)
        yield v


def check_case(N, r):
    u = min(r, 1 - r)
    eta = Q(1, 6) * (u / 72) ** N
    margin = (u / 4) ** N
    assert 0 < eta <= Q(1, 2 * N)
    assert 3 * eta * 18 ** N == margin / 2
    assert Q(factorial(N), N ** N) >= Q(1, 3 ** N)
    assert Q(1, 2 ** N) * Q(factorial(N), N ** N) >= Q(1, 6 ** N)
    v0 = reference(N, r)
    B0 = matrices(v0)
    J0 = jacobi(v0)
    minimum = Q(1)
    for m in range(N + 1):
        assert norm(B0[m]) == 1
        for n in range(N + 1):
            assert sum(B0[m][n]) == 1
            for k in range(N + 1):
                g = B0[m][n][k]
                assert g == tensor_coefficient(N, r, m, n, k)
                if triangle(m, n, k, N):
                    assert g >= margin
                    minimum = min(minimum, g)
                else:
                    assert g == 0
    # Fidelity to P1 Lemma A.1, independently obtained by two formulas above.
    assert B0[N][N][N] == (1 - r) ** N
    count = 0
    for v in samples(v0, eta):
        count += 1
        admissible(v)
        assert max(abs(x - y) for a, b in zip(v, v0) for x, y in zip(a, b)) <= eta
        J = jacobi(v)
        assert norm(J) == 1 and norm(sub(J, J0)) <= 3 * eta
        B = matrices(v)
        Delta = [sub(A, A0) for A, A0 in zip(B, B0)]
        D = [norm(A) for A in Delta]
        d = [Q(1)]
        for m in range(N):
            assert v[m][0] >= v0[m][0] / 2
            d.append(d[-1] * v[m][0])
        E = [x * y for x, y in zip(d, D)]
        assert E[0] == 0 and E[1] <= 3 * eta
        for m in range(1, N):
            a, b, c = v[m]
            da, db, dc = (v[m][i] - v0[m][i] for i in range(3))
            L = sub(left_mul(v, Delta[m]), scale(b, Delta[m]))
            L = sub(L, scale(c, Delta[m - 1]))
            # F is the exact inhomogeneous part of the subtracted recurrence.
            F = sub(left_mul(v, B0[m]), left_mul(v0, B0[m]))
            F = sub(sub(sub(F, scale(db, B0[m])), scale(dc, B0[m - 1])),
                    scale(da, B0[m + 1]))
            assert scale(a, Delta[m + 1]) == sub(L, scale(-1, F))
            assert norm(F) <= 6 * eta
            assert E[m + 1] <= 2 * E[m] + E[m - 1] + 6 * eta
        for m in range(N + 1):
            assert 0 < d[m] <= 1 and d[m] >= Q(1, 6 ** N)
            assert E[m] <= 3 ** (m + 1) * eta
            assert D[m] <= 3 * eta * 18 ** N
            for n in range(N + 1):
                assert sum(B[m][n]) == 1
                for k in range(N + 1):
                    g = B[m][n][k]
                    assert g >= margin / 2 if triangle(m, n, k, N) else g == 0
    print(f'N={N} r={r} arrays={count} eta={eta} reference_min={minimum} '
          f'coefficient_comparisons={count * (N + 1) ** 3} all_invariants=OK')
    return count


def main():
    envelope = [0, 3]
    for m in range(1, 8):
        envelope.append(2 * envelope[m] + envelope[m - 1] + 6)
    assert all(value <= 3 ** (m + 1) for m, value in enumerate(envelope))
    print('normalized_error_envelope=' + ','.join(map(str, envelope)))
    total = sum(check_case(N, r) for r in [Q(1, 2), Q(1, 3)] for N in range(1, 9))
    print(f'total_parameter_pairs=16 total_sample_arrays={total} '
          'fidelity_corner=OK tensor_formula=OK triangle_zeros=OK '
          'factorial_bound=OK scaled_recurrence=OK final_constants=OK')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
