#!/usr/bin/env python3
"""Exact checks for Zhao's fixed-dimension Recht-Re question; no numerical evidence."""
from itertools import combinations
from random import Random
import sympy as sp

R = sp.Rational


def clean(A):
    return A.applyfunc(sp.simplify)


def same(A, B):
    assert clean(A - B) == sp.zeros(A.rows, A.cols)


def total(B):
    return sum(B, sp.zeros(B[0].rows))


def ordered_q(B):
    return sum((A * C for i, A in enumerate(B)
                for j, C in enumerate(B) if i != j), sp.zeros(B[0].rows))


def q(B):
    T = total(B)
    return clean(T * T - sum((A * A for A in B), sp.zeros(T.rows)))


def c(k):
    assert k >= 1
    return R(k - 1, 4 * k)


def psd_rational(A):
    """Exact PSD criterion: symmetry and nonnegative principal minors of every order."""
    assert all(a.is_Rational for a in A)
    assert A == A.T
    for r in range(1, A.rows + 1):
        for inds in combinations(range(A.rows), r):
            assert A.extract(inds, inds).det() >= 0


def canonical(k):
    I = sp.eye(k)
    P = sp.ones(k) / k
    C = I - (1 - sp.sqrt(2) / 2) * P
    B = [clean(C[:, i] * C[:, i].T) for i in range(k)]
    return B, P, C


def embed(A, d):
    return sp.diag(A, sp.zeros(d - A.rows))


def fidelity():
    # Source Proposition 4.1: canonical Gram, sum, and ordered-product identities.
    for n in range(2, 10):
        B, P, C = canonical(n)
        I = sp.eye(n)
        same(P * P, P)
        same(C.T * C, I - P / 2)
        same(total(B), I - P / 2)
        for i, A in enumerate(B):
            same(A, A.T)
            same(A * A, (1 - R(1, 2 * n)) * A)
            for j in range(n):
                assert sp.simplify((C[:, i].T * C[:, j])[0]
                                   - (int(i == j) - R(1, 2 * n))) == 0
        expected = (I - P) / (2 * n) - c(n) * P
        same(q(B), expected)
        if n <= 4:
            same(ordered_q(B), expected)
        # Reproduce the source's A_i = n B_i scaling as well.
        same(q([n * A for A in B]), n * n * expected)
    print('FIDELITY_OK: source canonical identities and scaling for n=2..9; direct ordered sums for n=2..4')


def attainment():
    cases = 0
    for n in range(2, 10):
        for d in range(1, 9):
            k = min(n, d)
            base, P, C = canonical(k)
            B = [embed(A, d) for A in base] + [sp.zeros(d) for _ in range(n - k)]
            J, P = embed(sp.eye(k), d), embed(P, d)
            T, Q = clean(total(B)), q(B)
            same(T, J - P / 2)
            # Feasibility: each B is an outer product, and I-S=(I-J)+P/2.
            same(sp.eye(d) - T, (sp.eye(d) - J) + P / 2)
            same((sp.eye(d) - J)**2, sp.eye(d) - J)
            expected = (J - P) / (2 * k) - c(k) * P
            same(Q, expected)
            # Exact PSD certificate: Q+cI=(J-P)/(2k)+c(I-P).
            same(Q + c(k) * sp.eye(d), (J - P) / (2 * k) + c(k) * (sp.eye(d) - P))
            same((J - P)**2, J - P)
            same((sp.eye(d) - P)**2, sp.eye(d) - P)
            e = sp.Matrix([1 / sp.sqrt(k)] * k + [0] * (d - k))
            assert sp.simplify((e.T * e)[0]) == 1
            same(Q * e, -c(k) * e)
            if d < n:
                assert c(n) - c(d) == R(n - d, 4 * n * d) > 0
            # All-positive-definite approximants, checked at two exact epsilons.
            for eps in [R(1, 7), R(1, 31)]:
                Be = [clean((A + eps * sp.eye(d)) / (1 + n * eps)) for A in B]
                same(sp.eye(d) - total(Be), (sp.eye(d) - T) / (1 + n * eps))
                target = (-c(k) + (n - 1) * eps + n * (n - 1) * eps**2) / (1 + n * eps)**2
                same(q(Be) * e, target * e)
            cases += 1
    print(f'ATTAINMENT_OK: {cases} exact (n,d) pairs; includes k=1, d<n, d=n, d>n; 128 exact PD perturbations')


def merging_checks():
    rng = Random(261000155)
    cases = merges = dropped = 0
    for n in range(2, 7):
        for d in range(1, 5):
            for sample in range(3):
                vs = [sp.Matrix([rng.randrange(-2, 3) for _ in range(d)]) for _ in range(n)]
                # Every fourth family includes an actual annihilator.
                if sample == 0:
                    vs[-1] = sp.zeros(d, 1)
                denom = 1 + sum((v.T * v)[0] for v in vs)
                B = [v * v.T / denom for v in vs]
                T = total(B)
                for A in B:
                    psd_rational(A)
                psd_rational(sp.eye(d) - T)
                Q = q(B)
                same(Q, ordered_q(B))
                psd_rational(Q + c(min(n, d)) * sp.eye(d))
                x = sp.Matrix([rng.randrange(-2, 3) for _ in range(d)])
                if x == sp.zeros(d, 1):
                    x[0] = 1
                norm2 = (x.T * x)[0]
                before = (x.T * Q * x)[0]
                active = [A for A in B if A * x != sp.zeros(d, 1)]
                dropped += n - len(active)
                current = (x.T * q(active) * x)[0] if active else sp.Integer(0)
                assert current == before
                while len(active) > d:
                    z = [A * x for A in active]
                    assert all((x.T * v)[0] > 0 for v in z)
                    pairs = [(i, j) for i in range(len(active))
                             for j in range(i + 1, len(active)) if (z[i].T * z[j])[0] > 0]
                    assert pairs, 'obtuse-vector halfspace consequence failed'
                    i, j = pairs[0]
                    dot = (z[i].T * z[j])[0]
                    A = active[i] + active[j]
                    active = [M for h, M in enumerate(active) if h not in (i, j)] + [A]
                    after = (x.T * q(active) * x)[0]
                    assert after == current - 2 * dot < current
                    current = after
                    merges += 1
                assert len(active) <= min(n, d)
                assert before >= current >= -c(min(n, d)) * norm2
                if active:
                    psd_rational(sp.eye(d) - total(active))
                    r = len(active)
                    if r == 1:
                        assert current == 0
                    else:
                        assert current >= -c(r) * norm2
                    # Source Lemma 3.2: noncommutative lower-bound slack identity.
                    Tr = total(active)
                    D = sum((active[i] + active[j] - (active[i] - active[j])**2
                             for i in range(r) for j in range(i + 1, r)), sp.zeros(d))
                    same(q(active) + c(r) * sp.eye(d),
                         R(r - 1, r) * (Tr - sp.eye(d)/2)**2 + D/r)
                cases += 1
    assert merges > 0 and dropped > 0
    print(f'MERGING_OK: {cases} rational feasible families; {merges} strict merges; {dropped} annihilators removed; all principal-minor bounds exact')


def main():
    fidelity()
    attainment()
    merging_checks()
    print('SCOPE: finite exact checks and algebraic certificates; no proof of the universal dimension bound or Lean proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
