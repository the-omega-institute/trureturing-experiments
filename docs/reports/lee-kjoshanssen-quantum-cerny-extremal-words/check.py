#!/usr/bin/env python3
"""Exact finite checks for arXiv:2609.40154v1, section 6, item 4.

Python 3 and SymPy. No numerical evaluation, floating point, random search,
or assumption of the desired lower bound is used. Finite checks support the
preregistration; the all-length argument is in the draft, not proved by this run.
"""
from itertools import product
from math import isqrt
import sympy as S


def ceil_sqrt(n):
    r = isqrt(n)
    return r if r * r == n else r + 1


def words_at_most(m):
    for k in range(m + 1):
        yield from product((0, 1), repeat=k)


def transition_image(m, kind, word):
    """All columns of the word map; -1 denotes the zero vector."""
    image = list(range(m))
    for a in word:
        for j, q in enumerate(image):
            if q == -1:
                continue
            if kind == 'constant':
                image[j] = (-1 if q == 0 else q - 1) if a == 0 else q
            elif a == 0:
                image[j] = min(q + 1, m - 1)
            else:
                image[j] = -1 if q == m - 1 else 0
    return tuple(image)


def matrix_from_images(images):
    m = len(images)
    return S.Matrix(m, m, lambda i, j: int(images[j] == i))


def clean(X):
    return X.applyfunc(S.simplify)


def assert_zero(X):
    assert clean(X) == S.zeros(*X.shape), X


def basis(d):
    """Hilbert--Schmidt orthonormal traceless Hermitian matrices."""
    out = []
    for a in range(d):
        for b in range(a + 1, d):
            x, y = S.zeros(d), S.zeros(d)
            x[a, b] = x[b, a] = 1
            y[a, b], y[b, a] = -S.I, S.I
            out.extend((x / S.sqrt(2), y / S.sqrt(2)))
    for k in range(1, d):
        out.append(S.diag(*([1] * k + [-k] + [0] * (d - k - 1))) / S.sqrt(k * (k + 1)))
    assert len(out) == d * d - 1
    for i, g in enumerate(out):
        assert g.H == g and S.trace(g) == 0
        for j, h in enumerate(out):
            assert S.simplify(S.trace(g * h)) == int(i == j)
    return out


def choi(f, d):
    J = S.zeros(d * d)
    for a in range(d):
        for b in range(d):
            E = S.zeros(d)
            E[a, b] = 1
            J += S.kronecker_product(E, f(E))
    return clean(J)


def positive_ldl(A):
    """Exact LDL* certificate, strictly positive pivots, exact reconstruction."""
    assert_zero(A - A.H)
    n = A.rows
    L = S.eye(n)
    pivots = []
    for k in range(n):
        pivot = S.simplify(A[k, k] - sum(L[k, j] * S.conjugate(L[k, j]) * pivots[j] for j in range(k)))
        assert pivot.is_positive is True, pivot
        pivots.append(pivot)
        for i in range(k + 1, n):
            L[i, k] = S.simplify((A[i, k] - sum(L[i, j] * S.conjugate(L[k, j]) * pivots[j] for j in range(k))) / pivot)
    assert_zero(A - L * S.diag(*pivots) * L.H)
    return len(pivots)


def check_channels(m, kind):
    d = ceil_sqrt(m + 1)
    g = basis(d)[:m]
    eps = S.Rational(1, 2 * d * (m + 1))
    mixed = S.eye(d) / d
    start_index = m - 1 if kind == 'constant' else 0
    rho = mixed + eps * g[start_index]
    assert S.trace(rho) == 1
    positive_ldl(rho)
    channels = []
    for a in (0, 1):
        images = transition_image(m, kind, (a,))
        def T(X, images=images):
            return clean(sum((S.trace(g[j] * X) * g[q] for j, q in enumerate(images) if q != -1), S.zeros(d)))
        def phi(X, T=T):
            return clean(S.trace(X) * mixed + eps * T(X))
        JT, JP = choi(T, d), choi(phi, d)
        count = sum(q != -1 for q in images)
        assert S.simplify(S.trace(JT.H * JT)) == count
        assert_zero(JP - S.eye(d * d) / d - eps * JT)
        assert eps * eps * count <= S.Rational(1, 4 * d * d)
        # Stronger direct check: J(phi) >= I/(2d), certified by exact LDL*.
        positive_ldl(JP - S.eye(d * d) / (2 * d))
        for u in range(d):
            for v in range(d):
                assert S.simplify(sum(JP[u*d+k, v*d+k] for k in range(d))) == int(u == v)
        assert_zero(phi(mixed) - mixed)
        for j in range(m):
            target = S.zeros(d) if images[j] == -1 else eps * g[images[j]]
            assert_zero(phi(mixed + g[j]) - mixed - target)
        channels.append(phi)
    word = (0,) * m if kind == 'constant' else (0,) * (m - 1) + (1,)
    for a in word:
        rho = channels[a](rho)
    assert_zero(rho - mixed)
    # The actual reachable differences contain a nonzero multiple of every g_j:
    # the target mixed is reached by word, and its zero-prefix states are below.
    rho = mixed + eps * g[start_index]
    for j in range(m):
        idx = m - 1 - j if kind == 'constant' else j
        assert_zero(rho - mixed - eps ** (j + 1) * g[idx])
        rho = channels[0](rho)
    print(f"CPTP_EXACT kind={kind} m={m} d={d} eps={eps} choi_margin=1/{2*d} reachable_span={m}")


def main():
    total = 0
    for m in range(1, 13):
        target = (0,) * m
        zeros = [w for w in words_at_most(m) if all(q == -1 for q in transition_image(m, 'constant', w))]
        assert zeros == [target]
        # Source Theorem 3.1: nilpotency index exactly m, dimension obstruction m.
        J = matrix_from_images(transition_image(m, 'constant', (0,)))
        assert [int((J ** k).rank()) for k in range(m + 1)] == list(range(m, -1, -1))
        d = ceil_sqrt(m + 1)
        assert (d - 1) ** 2 - 1 < m <= d * d - 1
        print(f"FIDELITY constant m={m} unique_shortest={''.join(map(str,target))} minimal_dimension={d}")
        total += 2 ** (m + 1) - 1
    for m in range(2, 13):
        target = (0,) * (m - 1) + (1,)
        zeros = []
        for w in words_at_most(m):
            image = transition_image(m, 'marker', w)
            if all(q == -1 for q in image):
                zeros.append(w)
            if w != target:
                assert image[0] != -1
        assert zeros == [target]
        T = matrix_from_images(transition_image(m, 'marker', (0,)))
        B = matrix_from_images(transition_image(m, 'marker', (1,)))
        ranks = [int((T ** j).rank()) for j in range(m + 1)]
        assert ranks == list(range(m, 0, -1)) + [1]
        assert B * T ** (m - 1) == S.zeros(m)
        assert B * T ** (m - 2) != S.zeros(m)
        d = ceil_sqrt(m + 1)
        assert (d - 1) ** 2 - 1 < m <= d * d - 1
        print(f"FAMILY m={m} unique_shortest={''.join(map(str,target))} image_ranks={ranks} dimension_threshold={d}")
        total += 2 ** (m + 1) - 1
    for m in range(2, 9):
        check_channels(m, 'marker')
    check_channels(8, 'constant')
    print(f"SOURCE_VALUE qc(0^8)=3; FAMILY_VALUE predicted qc(0001)=3; enumerated_word_tests={total}")
    print("SCOPE finite exact checks; universal lower bound and all-length settlement await Lean proof")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
