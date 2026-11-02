from fractions import Fraction as Q
from math import isqrt


def qsqrt(x: Q) -> Q:
    if x < 0:
        raise AssertionError("negative square root")
    rn, rd = isqrt(x.numerator), isqrt(x.denominator)
    if rn * rn != x.numerator or rd * rd != x.denominator:
        raise AssertionError(f"non-square rational {x}")
    return Q(rn, rd)


def power(x: Q, alpha: Q) -> Q:
    if alpha == 0:
        return Q(1)
    if alpha == 1:
        return x
    if alpha == Q(1, 2):
        return qsqrt(x)
    raise AssertionError(f"unsupported exact exponent {alpha}")


def cadd(z, w):
    return (z[0] + w[0], z[1] + w[1])


def csub(z, w):
    return (z[0] - w[0], z[1] - w[1])


def cscale(x, z):
    return (x * z[0], x * z[1])


def norm_sq(z):
    return z[0] * z[0] + z[1] * z[1]


def spectral_quantities(lambdas, alpha, matrix):
    n = len(lambdas)
    mean = (Q(0), Q(0))
    for i in range(n):
        mean = cadd(mean, cscale(lambdas[i], matrix[i][i]))
    centered = [[csub(matrix[i][j], mean if i == j else (Q(0), Q(0)))
                 for j in range(n)] for i in range(n)]
    weights = [[norm_sq(centered[i][j]) for j in range(n)] for i in range(n)]
    pa = [power(x, alpha) for x in lambdas]
    pb = [power(x, Q(1) - alpha) for x in lambdas]
    v = Q(0)
    info = Q(0)
    joint = Q(0)
    for i in range(n):
        for j in range(n):
            w = weights[i][j]
            s = lambdas[i] + lambdas[j]
            cross = pa[i] * pb[j] + pa[j] * pb[i]
            plus = (pa[i] + pa[j]) * (pb[i] + pb[j])
            v += Q(1, 2) * s * w
            info += Q(1, 2) * (s - cross) * w
            joint += Q(1, 2) * plus * w
    return v, info, joint


def check_case(lambdas, alpha, matrix):
    v, info, joint = spectral_quantities(lambdas, alpha, matrix)
    m, M = min(lambdas), max(lambdas)
    S, T = M + m, power(M, alpha) * power(m, Q(1) - alpha) + power(M, Q(1) - alpha) * power(m, alpha)
    D, P = S - T, S + T
    c = D / S
    assert joint == 2 * v - info
    assert Q(0) <= info <= c * v
    assert c * v - info >= 0
    assert (2 - c) * v - info >= 0
    # The factorization is the new scalar comparison after J = 2V - I.
    lhs_gap = c * c * v * v - (c / (2 - c)) * info * joint
    rhs_gap = (c / (2 - c)) * (c * v - info) * ((2 - c) * v - info)
    assert lhs_gap == rhs_gap
    assert lhs_gap >= 0
    # Avoid square roots as evidence: all factors are nonnegative, so this is
    # exactly equivalent to sqrt(D/P) U <= D/(M+m) V.
    assert c * c * v * v >= (D / P) * info * joint
    return v, info, joint, D, P


def main():
    # Fidelity to the source's explicit two-dimensional diagonal calculation.
    lambdas = [Q(9, 25), Q(16, 25)]
    matrix = [[(Q(1, 3), Q(1, 5)), (Q(2, 3), Q(-1, 4))],
              [(Q(-3, 5), Q(2, 7)), (Q(4, 9), Q(-1, 6))]]
    v, info, joint, D, P = check_case(lambdas, Q(1, 2), matrix)
    a, b = matrix[0][0], matrix[1][1]
    cxy, dxy = matrix[0][1], matrix[1][0]
    mean_gap = norm_sq(csub(a, b))
    source_v = Q(1, 2) * (norm_sq(cxy) + norm_sq(dxy)) + lambdas[0] * lambdas[1] * mean_gap
    assert v == source_v
    print("fidelity: two-dimensional NV formula and comparison reproduced")

    cases = [
        ([Q(9, 25), Q(16, 25)], Q(1, 2)),
        ([Q(1, 9), Q(4, 9), Q(4, 9)], Q(1, 2)),
        ([Q(1, 25), Q(4, 25), Q(4, 25), Q(16, 25)], Q(1, 2)),
        ([Q(9, 25), Q(16, 25)], Q(0)),
        ([Q(9, 25), Q(16, 25)], Q(1)),
    ]
    matrices = [
        matrix,
        [[(Q(1, 7), Q(1, 5)), (Q(2, 7), Q(0)), (Q(-1, 3), Q(1, 8))],
         [(Q(3, 11), Q(-1, 9)), (Q(2, 5), Q(1, 6)), (Q(1, 4), Q(-2, 7))],
         [(Q(-2, 9), Q(1, 10)), (Q(1, 6), Q(-1, 11)), (Q(3, 8), Q(1, 12))]],
        [[(Q(1, 7), Q(1, 5)), (Q(2, 7), Q(0)), (Q(-1, 3), Q(1, 8)), (Q(1, 9), Q(1, 13))],
         [(Q(3, 11), Q(-1, 9)), (Q(2, 5), Q(1, 6)), (Q(1, 4), Q(-2, 7)), (Q(-2, 9), Q(1, 10))],
         [(Q(-2, 9), Q(1, 10)), (Q(1, 6), Q(-1, 11)), (Q(3, 8), Q(1, 12)), (Q(2, 13), Q(-1, 14))],
         [(Q(1, 5), Q(-1, 15)), (Q(-1, 7), Q(1, 16)), (Q(2, 11), Q(1, 17)), (Q(1, 3), Q(-1, 18))]],
        matrix,
        matrix,
    ]
    for (lambdas, alpha), mat in zip(cases, matrices):
        check_case(lambdas, alpha, mat)
    print(f"exact checks: {len(cases)} spectra/exponents, including n=2,3,4")

    # The source's sharp endpoint family A = E_(1,n) attains equality.
    for lambdas in cases[:3] and [cases[0][0], cases[1][0], cases[2][0]]:
        n = len(lambdas)
        e = [[(Q(0), Q(0)) for _ in range(n)] for _ in range(n)]
        e[0][n - 1] = (Q(1), Q(0))
        v, info, joint, D, P = check_case(lambdas, Q(1, 2), e)
        S = max(lambdas) + min(lambdas)
        assert v == S / 2 and info == D / 2 and joint == P / 2
        assert (D / P) * info * joint == (D / S) ** 2 * v * v
    print("sharpness: E_(1,n) equality verified for n=2,3,4")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
