#!/usr/bin/env python3
"""Exact finite checks for arXiv:2611.04009 Appendix B; not a proof of its a.e. claim.

All certificates use rational algebra and rational upper/lower bounds for ln.
Random complex rational Cayley unitaries are test cases, not Haar samples.
An optional independent float diagnostic is provided separately.
"""
from fractions import Fraction as F
from functools import lru_cache
from random import Random
import sympy as s


def rat(x):
    x = s.cancel(x)
    assert x.is_Rational, x
    return F(int(x.p), int(x.q))


def sq(x):
    return rat(s.expand(x * s.conjugate(x)))


@lru_cache(None)
def log_bounds(x):
    """ln(x) enclosed by 2 atanh((x-1)/(x+1)), with an exact geometric tail."""
    x = F(x)
    assert x > 0
    k = 0
    while x >= 2:
        x /= 2
        k += 1
    while x < 1:
        x *= 2
        k -= 1

    def base(y):
        z = (y - 1) / (y + 1)
        z2 = z * z
        term = z
        total = F(0)
        for n in range(48):
            total += 2 * term / (2 * n + 1)
            term *= z2
        tail = 2 * term / (97 * (1 - z2))
        return total, total + tail

    lo, hi = base(x)
    l2, u2 = base(F(2))
    if k >= 0:
        return lo + k * l2, hi + k * u2
    return lo + k * u2, hi + k * l2


def entropy_bounds(p):
    assert all(x >= 0 for x in p) and sum(p) == 1
    lo = hi = F(0)
    for x in p:
        if x:
            a, b = log_bounds(x)
            lo -= x * b
            hi -= x * a
    return lo, hi


def interval_mul(a, b):
    p = [x * y for x in a for y in b]
    return min(p), max(p)


def variance_bounds(p):
    logs = [log_bounds(x) for x in p]
    m = (sum(x * z[0] for x, z in zip(p, logs)),
         sum(x * z[1] for x, z in zip(p, logs)))
    m2 = interval_mul(m, m)
    second = [interval_mul(z, z) for z in logs]
    return (sum(x * z[0] for x, z in zip(p, second)) - m2[1],
            sum(x * z[1] for x, z in zip(p, second)) - m2[0])


def kd(U, rho):
    return s.Matrix(U.rows, U.cols,
                    lambda i, j: s.simplify(s.conjugate(U[i, j]) * (rho * U)[i, j]))


def fidelity():
    U = s.Matrix([[-1, 2, 2], [2, -1, 2], [2, 2, -1]]) / 3
    assert U.H * U == s.eye(3)
    v = s.Matrix([1, -1, 0]) / s.sqrt(2)
    Q = kd(U, v * v.H)
    assert Q == s.Matrix([[1, 2, 0], [2, 1, 0], [0, 0, 0]]) / 6
    assert [sum(Q.row(i)) for i in range(3)] == [s.Rational(1, 2), s.Rational(1, 2), 0]
    assert [sum(Q.col(j)) for j in range(3)] == [s.Rational(1, 2), s.Rational(1, 2), 0]
    for i in range(3):
        a = s.eye(3)[:, i]
        Qa = kd(U, a * a.H)
        row = [sq(U[i, j]) for j in range(3)]
        assert sorted(row) == [F(1, 9), F(4, 9), F(4, 9)]
        assert [rat(Qa[i, j]) for j in range(3)] == row
        assert all(Qa[k, j] == 0 for k in range(3) if k != i for j in range(3))
        assert entropy_bounds([rat(z) for z in Qa]) == entropy_bounds(row)
    U0 = s.Matrix([[-8, 15, -6, -6], [3, -8, -12, -12],
                   [-12, -6, 10, -9], [-12, -6, -9, 10]]) / 19
    assert U0.H * U0 == s.eye(4)
    psi = s.Matrix([2, -1, 0, 0]) / s.sqrt(5)
    assert U0.H * psi == s.Matrix([-1, 2, 0, 0]) / s.sqrt(5)
    assert kd(U0, psi * psi.H) == s.Matrix([[16, 60, 0, 0], [3, 16, 0, 0],
                                                         [0, 0, 0, 0], [0, 0, 0, 0]]) / 95
    print('FIDELITY_OK: Appendix B matrices (58) and U_0, all displayed KD coordinates and marginals')


def random_unitary(d, rng):
    X = s.Matrix(d, d, lambda i, j: s.Rational(rng.randint(-3, 3), 7)
                 + s.I * s.Rational(rng.randint(-3, 3), 7))
    K = X - X.H
    U = (s.eye(d) - K) * (s.eye(d) + K).inv()
    U = U.applyfunc(s.cancel)
    assert (U.H * U).applyfunc(s.simplify) == s.eye(d)
    return U


def kd_real_nullity(U):
    d = U.rows
    basis = []
    for i in range(d):
        H = s.zeros(d); H[i, i] = 1; basis.append(H)
    for i in range(d):
        for j in range(i + 1, d):
            H = s.zeros(d); H[i, j] = H[j, i] = 1; basis.append(H)
            H = s.zeros(d); H[i, j] = s.I; H[j, i] = -s.I; basis.append(H)
    columns = [s.Matrix([s.im(z) for z in kd(U, H)]) for H in basis]
    return d * d - s.Matrix.hstack(*columns).rank()


def check_case(U, case):
    d = U.rows
    assert all(z != 0 for z in U)
    assert kd_real_nullity(U) == 2 * d - 1
    rows = [[sq(U[i, j]) for j in range(d)] for i in range(d)]
    cols = [[sq(U[i, j]) for i in range(d)] for j in range(d)]
    laws = rows + cols
    bounds = [entropy_bounds(p) for p in laws]
    vtx = min(range(2 * d), key=lambda k: bounds[k][1])
    for k in range(2 * d):
        assert sorted(laws[k]) == sorted(laws[vtx]) or bounds[vtx][1] < bounds[k][0]
    if vtx >= d:
        U = U.H
        i = vtx - d
    else:
        i = vtx
    p = [sq(U[i, j]) for j in range(d)]
    assert len(set(p)) > 1
    h = entropy_bounds(p)
    V = variance_bounds(p)
    assert V[0] > 0
    # Exact symbolic identity for the actual logarithmic tangent, no numerical logs.
    z = s.symbols('z:' + str(d), real=True)
    pp = [s.Rational(x.numerator, x.denominator) for x in p]
    m = sum(x * y for x, y in zip(pp, z))
    var = sum(x * (y - m)**2 for x, y in zip(pp, z))
    derivative = -sum(2 * x * (y - m) * (1 + y) for x, y in zip(pp, z))
    assert s.expand(derivative + 2 * var) == 0
    # A rational approximation to L, accurate within 2^-63 in each diagonal entry.
    scale = 2**64
    L = []
    for x in p:
        lo, hi = log_bounds(x)
        assert hi - lo < F(1, scale)
        L.append(s.Rational(((lo + hi) * scale / 2).__floor__(), scale))
    a = s.eye(d)[:, i]
    m0 = sum(s.Rational(x.numerator, x.denominator) * y for x, y in zip(p, L))
    v = (U * s.diag(*L) * U.H * a - m0 * a).applyfunc(s.simplify)
    assert (a.H * v)[0] == 0
    Vr = sum(sq(w) for w in v)
    assert Vr > 0
    def entropies(t):
        x = a + s.Rational(t.numerator, t.denominator) * v
        norm = sum(sq(w) for w in x)
        assert norm == 1 + t*t*Vr
        pa = [sq(w) / norm for w in x]
        pb = [sq(w) / norm for w in U.H * x]
        return entropy_bounds(pa), entropy_bounds(pb)
    # The rational state provides a certified strict witness below the global vertex minimum.
    t = F(1, 2**20)
    HA, HB = entropies(t)
    Fhi = HA[1] + HB[1]
    assert Fhi < h[0]
    # Finite-difference agreement with -2 Var_p(ln p), using rigorous rational intervals.
    diff = ((HB[0] - h[1]) / t, (HB[1] - h[0]) / t)
    residual = (diff[0] + 2 * V[0], diff[1] + 2 * V[1])
    assert max(abs(residual[0]), abs(residual[1])) < F(1, 1000)
    assert HA[1] / t < F(1, 1000)
    assert Fhi < h[0] - F(1, 2**60)
    print(f'EXACT_CASE_OK d={d} case={case} KD-real-nullity={2*d-1} '
          f'vertex={vtx} t=1/1048576 gap>2^-60 '
          'abs(B-finite-difference+2V)<1/1000 A(t)/t<1/1000')


def main():
    fidelity()
    rng = Random(261104009)
    for d in (2, 3, 4):
        for case in range(3):
            check_case(random_unitary(d, rng), case)
    print('All certificates above use exact arithmetic; finite cases do not establish Haar-a.e. quantifiers.')
    import sys
    if '--haar-diagnostics' in sys.argv:
        import haar_diagnostics
        haar_diagnostics.main()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
