#!/usr/bin/env python3
"""Exact checks of the n=3 Camargo--Nishida Page/ME formula.

Only integer and Fraction arithmetic is used. Determinants are assembled from
replica permutations, independently of the proposed nine-character spectrum.
No logarithm or floating-point approximation is evaluated.
"""
from fractions import Fraction as F
from itertools import product


def determinant(rows):
    a = [list(map(F, row)) for row in rows]
    value = F(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            value = -value
        p = a[j][j]
        value *= p
        for i in range(j + 1, len(a)):
            if a[i][j]:
                c = a[i][j] / p
                for k in range(j + 1, len(a)):
                    a[i][k] -= c * a[j][k]
                a[i][j] = F(0)
    return value


def replicas(tripartite):
    return list(product(range(3), repeat=2)) if tripartite else list(range(3))


def twist(party, replica, tripartite):
    if tripartite:
        i, j = replica
        return ((i + 1) % 3, j) if party == 0 else (
            (i, (j + 1) % 3) if party == 1 else replica)
    return (replica + 1) % 3 if party == 0 else replica


def replica_det(parts, r, tripartite):
    # Existing integral-to-replica-determinant reduction, no character formula.
    rs = replicas(tripartite)
    n = sum(parts)
    q = [[sum(F(size, n) for party, size in enumerate(parts)
              if row == twist(party, col, tripartite))
          for col in rs] for row in rs]
    m = len(rs)
    gram = [[sum(q[k][i] * q[k][j] for k in range(m))
             for j in range(m)] for i in range(m)]
    h = [[F(i == j) * ((1 + r) / 2) ** 2
          - ((1 - r) / 2) ** 2 * gram[i][j]
          for j in range(m)] for i in range(m)]
    return determinant(h)


def full_det(parts, r):
    # Direct 9N x 9N Gaussian precision after the scale t is cancelled.
    rs = replicas(True)
    modes = [party for party, size in enumerate(parts) for _ in range(size)]
    indices = [(rep, k) for rep in rs for k in range(len(modes))]
    n = len(modes)
    return determinant([
        [F(i == j) - (1 - r) / (2 * n) * (
            int(i[0] == j[0]) + int(
                twist(modes[i[1]], i[0], True) ==
                twist(modes[j[1]], j[0], True)))
         for j in indices] for i in indices])


def source_parameters(n, r):
    # t^2 is rational; a=t*(N-1+r)/N and e=t*(r-1)/N.
    t2 = (n - 1 + 1 / r) / (n - 1 + r)
    a_over_t = (n - 1 + r) / n
    e_over_t = (r - 1) / n
    a2 = t2 * a_over_t ** 2
    d = a2 - 1
    assert t2 > 0 and a_over_t > 0 and d >= 0
    # Exact sign and square identify the source's nonnegative radical product.
    radical = F(n - 1, n) * (1 - r * r) / r
    assert radical >= 0
    assert radical ** 2 == d * (d * n * n + 4 * (n - 1))
    assert (d * (n - 2) - radical) / (
        2 * a_over_t * t2 * (n - 1)) == e_over_t
    assert (1 - r) ** 2 / r == F(n * n, n - 1) * d
    return a2, d


def proposed_ratio(parts, d):
    a, b, c = parts
    u = 3 * d / (4 * (a + b + c - 1))
    return (1 + u * (a * b + a * c + b * c)) ** 2 / (
        (1 + u * a * (b + c)) * (1 + u * b * (a + c)) *
        (1 + u * c * (a + b)))


def main():
    fidelity = small = full = 0
    r_values = [F(1), F(1, 2), F(1, 3), F(1, 5)]
    for n in range(3, 10):
        for a in range(1, n - 1):
            for b in range(1, n - a):
                c = n - a - b
                parts = (a, b, c)
                for r in r_values:
                    a2, d = source_parameters(n, r)
                    dt = replica_det(parts, r, True)
                    db = [replica_det((p, n - p), r, False) for p in parts]
                    assert dt > 0 and all(x > 0 for x in db)
                    ratio = proposed_ratio(parts, d)
                    assert ratio > 0
                    # exp(24 GM) = r^9 det(H_trip)^2 / prod(det(H_bip))^3.
                    # Squaring is harmless: both exp(12 GM) and ratio are positive.
                    assert r ** 9 * dt ** 2 == ratio ** 2 * (db[0] * db[1] * db[2]) ** 3
                    x, y, z = [F(p, n) for p in parts]
                    ds = [r + 3 * (1-r)**2 * p / 4 for p in
                          (x*(1-x), y*(1-y), z*(1-z), x*y+x*z+y*z)]
                    assert dt == r * (ds[0]*ds[1]*ds[2]*ds[3])**2
                    small += 1
                    # Three n=3 entries stated in the source's Table 1.
                    table = {
                        (1, 1, 1): (9*a2-1)**2/(3*a2+1)**3,
                        (2, 1, 1): (5*a2-1)**2/(a2*(3*a2+1)**2),
                        (3, 1, 1): (21*a2-5)**2/(2*(3*a2+1)**2*(9*a2-1)),
                    }
                    if parts in table:
                        assert ratio == table[parts]
                        fidelity += 1
                    if r == F(1, 2) and n <= 5:
                        assert full_det(parts, r) == dt
                        full += 1
    print(f'SOURCE_FIDELITY_OK: Table 1, three n=3 rows, {fidelity} exact instances')
    print(f'ALL_SMALL_CASES_OK: every ordered positive tripartition with 3<=N<=9, r in {r_values}; {small} cases')
    print(f'FULL_GAUSSIAN_DETERMINANTS_OK: every ordered positive tripartition with 3<=N<=5 at r=1/2; {full} cases')
    configs = 0
    for d in [F(1, 10), F(1), F(3), F(15)]:
        for k in range(1, 21):
            n = 4*k
            zp = 3*d*n*n/(64*(n-1))
            fp = (1+5*zp)**2/((1+4*zp)*(1+3*zp)**2)
            assert proposed_ratio((2*k,k,k), d) == fp
            assert fp*n*d/F(400,27) == F(n-1,n)*(1+1/(5*zp))**2/(
                (1+1/(4*zp))*(1+1/(3*zp))**2)
            n = 3*k
            zm = d*n*n/(12*(n-1))
            fm = (1+3*zm)**2/(1+2*zm)**3
            assert proposed_ratio((k,k,k), d) == fm
            assert fm*n*d/F(27,2) == F(n-1,n)*(1+1/(3*zm))**2/(1+1/(2*zm))**3
            configs += 2
    # Constants and absolute-error coefficient bookkeeping, without logs.
    assert F(25,36)*F(64,3) == F(400,27)
    assert F(9,8)*12 == F(27,2)
    assert F(400,27)/F(27,2) == F(800,729)
    assert 2*F(1,5)+F(1,4)+2*F(1,3) == F(79,60)
    assert 2*F(1,3)+3*F(1,2) == F(13,6)
    assert all(proposed_ratio(parts, F(0)) == 1
               for parts in [(1,1,1), (2,1,1), (2,2,2), (3,4,5)])
    print(f'PAGE_ME_FACTORIZATION_OK: {configs} exact configurations; asymptotic constants and error coefficients checked')
    print('BOUNDARY_OK: a=1 gives exp(12 GM)=1')
    print('LIMITATION: finite exact checks and algebraic factorization are not a proof of the quantified claim or logarithm bounds')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
