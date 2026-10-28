#!/usr/bin/env python3
"""Exact finite checks for Mandal--Mehatari--Das, arXiv:2205.00310v1.

Requires Python 3 and SymPy. Only integer/rational polynomial arithmetic and
certified rational real-root isolation are used. No numerical eigenvalues.
The finite ranges below do not prove the unbounded conjecture.
"""
from functools import lru_cache
from itertools import combinations
import sympy as sp

z = sp.Symbol("z")
EPS = sp.Rational(1, 10**16)


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def compositions(total, length):
    for cuts in combinations(range(1, total), length - 1):
        points = (0,) + cuts + (total,)
        yield tuple(points[i + 1] - points[i] for i in range(length))


def seidel(runs):
    word = tuple(bit for j, size in enumerate(runs) for bit in [j % 2] * size)
    n = len(word)
    return sp.Matrix(n, n, lambda i, j: 0 if i == j else
                     (-1 if word[min(i, j)] == 0 and word[max(i, j)] == 1 else 1))


def core(k):
    d = 2 * k - 1
    return sp.Matrix(d, d, lambda i, j:
                     -1 if min(i, j) % 2 == 0 and max(i, j) % 2 == 1 else 1)


def compressed_data(runs):
    k = len(runs) // 2
    sizes = (runs[0] + runs[-1],) + runs[1:-1]
    d, n = len(sizes), sum(runs)
    # First and last runs are merged after switching the last run.
    labels = [j if j < 2 * k - 1 else 0
              for j, size in enumerate(runs) for _ in range(size)]
    Z = sp.Matrix(n, d, lambda i, j: int(labels[i] == j))
    signs = [1] * (n - runs[-1]) + [-1] * runs[-1]
    P = sp.diag(*signs)
    D = sp.diag(*sizes)
    return core(k), D, Z, P


def polynomial(A):
    return sp.Poly(A.charpoly(z).as_expr(), z, domain=sp.QQ)


@lru_cache(maxsize=None)
def roots(coeffs):
    p = sp.Poly.from_list(list(coeffs), z, domain=sp.QQ)
    intervals = tuple(p.intervals(eps=EPS))
    require(sum(mult for _, mult in intervals) == p.degree(), "all roots real")
    return intervals


def root_intervals(p):
    return roots(tuple(p.all_coeffs()))


def energy_bounds(p):
    low = high = sp.Rational(0)
    for (a, b), mult in root_intervals(p):
        if a >= 0:
            low += mult * a
            high += mult * b
        elif b <= 0:
            low -= mult * b
            high -= mult * a
        else:
            raise AssertionError("root interval straddles zero")
    return low, high


def in_gamma(runs):
    return all(size == 1 for size in runs[1:-1])


def check_source_fidelity():
    # Lemma 4.12 states Spec_S(H) for H in Gamma_{n,2}. The exact
    # polynomial below encodes its eigenvalues and their multiplicities.
    count = 0
    for n in range(4, 15):
        expected = sp.Poly((z + 1)**(n - 3) * (z - 1) *
                           (z*z - (n - 4)*z - 3*n + 7), z, domain=sp.QQ)
        radical = sp.sqrt(n*n + 4*n - 12)
        positive = (n - 4 + radical) / 2
        negative = (n - 4 - radical) / 2
        require(positive > 0 and negative < 0, "Lemma 4.12 signs")
        # Therefore its energy is n-2 + sqrt(n^2+4n-12), as in Theorem 4.13.
        require(sp.simplify((n - 3) + 1 + positive - negative -
                            (n - 2 + radical)) == 0, "source minimum energy")
        for first in range(1, n - 2):
            runs = (first, 1, 1, n - 2 - first)
            require(polynomial(seidel(runs)) == expected, ("Lemma 4.12", runs))
            count += 1
    print(f"FIDELITY Lemma 4.12 and its energy consequence: {count} Gamma_{{n,2}} graphs, 4 <= n <= 14")


def check_inverse_and_rotation():
    count = 0
    for k in range(2, 13):
        M = core(k)
        d = M.rows
        N = sp.zeros(d)
        for i in range(d - 1):
            N[i, i + 1] = N[i + 1, i] = sp.Rational((-1)**(i + 1), 2)
        N[0, d - 1] = N[d - 1, 0] = sp.Rational(1, 2)
        require(M*N == sp.eye(d) and N*M == sp.eye(d), ("inverse", k))
        require(all(sum(abs(N[i, j]) for j in range(d)) == 1
                    for i in range(d)), ("absolute row sums", k))
        signs = [-1 if i == d - 2 else 1 for i in range(d)]
        for i in range(d):
            for j in range(d):
                require(M[(i + 2) % d, (j + 2) % d] == signs[i]*signs[j]*M[i, j],
                        ("signed rotation", k, i, j))
        count += 1
    print(f"IDENTITIES inverse signed cycle and signed two-step rotation: {count} cores, 2 <= k <= 12")


def check_all_small_cases():
    total = strict = equality = compressions = 0
    for n in range(2, 11):
        for k in range(1, n // 2 + 1):
            baseline = (n - 2*k + 1,) + (1,)*(2*k - 1)
            basepoly = polynomial(seidel(baseline))
            _, basehigh = energy_bounds(basepoly)
            for runs in compositions(n, 2*k):
                S = seidel(runs)
                p = polynomial(S)
                if in_gamma(runs):
                    # Source Lemma 4.10, also certifies exact energy equality.
                    require(p == basepoly, ("Gamma cospectrality", runs))
                    equality += 1
                else:
                    low, _ = energy_bounds(p)
                    require(low > basehigh, ("strict conjectured inequality", runs))
                    strict += 1
                M, D, Z, P = compressed_data(runs)
                d = M.rows
                require(P*(S + sp.eye(n))*P == Z*M*Z.T, ("compression", runs))
                require(Z.T*Z == D, ("class sizes", runs))
                q = polynomial(M*D)
                shifted = sp.Poly((z + 1)**(n - d) * q.as_expr().subs(z, z + 1), z, domain=sp.QQ)
                require(p == shifted, ("characteristic compression identity", runs))
                if k >= 2:
                    ri = root_intervals(q)
                    positive = sum(mult for (a, b), mult in ri if a >= 1)
                    negative = sum(mult for (a, b), mult in ri if b <= -1)
                    require(positive == k and negative == k - 1,
                            ("inertia and exclusion (-1,1)", runs))
                    # Exact charpoly shift + these certified signs imply
                    # SE(G) = n-2k + tr|D^(1/2) M D^(1/2)| in this case.
                    compressions += 1
                total += 1
    require(total == 511, "complete small-case enumeration")
    print(f"SMALL CASES all positive compositions of even length, 2 <= n <= 10: {total}")
    print(f"CONJECTURE exact energy comparison: strict={strict}, equality={equality}")
    print(f"COMPRESSION exact matrix/characteristic identities: {total}; energy-identity sign certificates (k>=2): {compressions}")
    print("ARITHMETIC integers, exact radicals and rational real-root isolation; no floating-point evidence")


def main():
    print("Mandal--Mehatari--Das, arXiv:2205.00310v1, Conjecture 4.14")
    print(f"SymPy {sp.__version__}; root-isolation width <= 1/10^16")
    check_source_fidelity()
    check_inverse_and_rotation()
    check_all_small_cases()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
