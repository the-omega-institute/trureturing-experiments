#!/usr/bin/env python3
"""Exact finite checks for the proposed full-matrix CPTP coproduct classification.

These checks validate finite families and proof identities, not the universal theorem.
All matrix entries belong to Q(i); there is no floating-point computation.
"""
import sympy as S

Q = S.Rational
I = S.I
checks = 0


def require(ok, message):
    global checks
    if not ok:
        raise AssertionError(message)
    checks += 1


def eq(a, b, message):
    require(a == b, message)


def zero(n, m=None):
    return S.zeros(n, n if m is None else m)


def unit(n, m, i, j):
    a = zero(n, m)
    a[i, j] = 1
    return a


def apply(ks, a, out):
    return sum((k * a * k.H for k in ks), zero(out))


def complete(ks, inp):
    return sum((k.H * k for k in ks), zero(inp))


def support(rho):
    cols = rho.columnspace()
    if not cols:
        return zero(rho.rows)
    b = S.Matrix.hstack(*cols)
    return b * (b.H * b).inv() * b.H


def rotation(c, a, b, imag=False):
    v = S.eye(c)
    v[a, a] = v[b, b] = Q(3, 5)
    if imag:
        v[a, b] = v[b, a] = Q(4, 5) * I
    else:
        v[a, b], v[b, a] = -Q(4, 5), Q(4, 5)
    eq(v.H * v, S.eye(c), 'rational Gaussian rotation unitary')
    return v


def channel_pair(n, m, extra, pure):
    r, s = (1, 1) if pure else (n, m)
    c = r + s + extra
    v = rotation(c, 0, r, imag=True)
    if c > 2:
        v = rotation(c, 1, c - 1) * v
    eq(v.H * v, S.eye(c), 'combined unitary')
    families = []
    for inp, offset in [(n, 0), (m, r)]:
        if pure:
            ks = [v * unit(c, inp, offset, j) for j in range(inp)]
        else:
            j = zero(c, inp)
            for t in range(inp):
                j[offset + t, t] = 1
            shift = zero(inp)
            for t in range(inp):
                shift[(t + 1) % inp, t] = 1
            ks = [Q(3, 5) * v * j, Q(4, 5) * v * j * shift]
        eq(complete(ks, inp), S.eye(inp), 'Kraus completeness for injection family')
        families.append(ks)
    return c, r, v, families


def fidelity():
    # The source uses dim_R CPTP(n,z) = n^2 (z^2 - 1), for n,z > 0.
    # Compute the real rank of TP constraints on Hermitian-to-Hermitian maps.
    # The depolarizing Choi matrix I_(nz)/z is strictly positive, so this
    # affine dimension is also the dimension of the CPTP convex body.
    count = 0
    for n in range(1, 5):
        for z in range(1, 5):
            traces = [1] * z + [0] * (z * z - z)
            constraints = S.zeros(n * n, n * n * z * z)
            for a in range(n * n):
                for b, tr in enumerate(traces):
                    constraints[a, a * z * z + b] = tr
            dimension = constraints.cols - constraints.rank()
            eq(dimension, n * n * (z * z - 1), 'source hom-set dimension')
            depolarizing_choi = Q(1, z) * S.eye(n * z)
            require(all(depolarizing_choi[t, t] > 0 for t in range(n * z)),
                    'strictly positive depolarizing Choi matrix')
            count += 1
    print(f'FIDELITY: source hom-set dimension formula, {count} positive (n,z) pairs, n,z <= 4')
    for n in range(1, 129):
        require(n * n < n * n + 1 < (n + 1) ** 2,
                'source Corollary 7.3 integer-square exclusion')
    print('FIDELITY: Corollary 7.3 arithmetic exclusion, 1 <= n <= 128')
    triples = [(n, m, c) for n in range(1, 14) for m in range(n, 14)
               for c in range(1, 20) if c * c == n * n + m * m]
    require((3, 4, 5) in triples and (5, 12, 13) in triples,
            'source Theorem 7.2 leaves Pythagorean pairs')
    print(f'FIDELITY: remaining square-dimension pairs up to 13: {triples}')


def obstruction(n, m, extra, pure):
    c, split, v, families = channel_pair(n, m, extra, pure)
    rhos, qs = [], []
    for inp, ks in zip((n, m), families):
        rho = apply(ks, Q(1, inp) * S.eye(inp), c)
        eq(S.trace(rho), 1, 'normalized image trace')
        require(rho != zero(c), 'normalized image nonzero')
        q = support(rho)
        eq(q.H, q, 'support projection Hermitian')
        eq(q * q, q, 'support projection idempotent')
        eq(q * rho, rho, 'support contains normalized output')
        for k in ks:
            eq(q * k, k, 'EVERY Kraus range lies in normalized-output support')
        qs.append(q)
        rhos.append(rho)
    eq(qs[0] * qs[1], zero(c), 'orthogonal supports')
    # The discriminator has one row of v^* per Kraus operator.
    ls = []
    for a in range(c):
        k = zero(2, c)
        k[0 if a < split else 1, :] = v[:, a].H
        ls.append(k)
    eq(complete(ls, c), S.eye(c), 'discriminator Kraus completeness')
    p0, p1 = S.diag(1, 0), S.diag(0, 1)
    e = sum((k.H * p0 * k for k in ls), zero(c))
    eq(e, qs[0], 'Heisenberg effect equals first support in this family')
    for j, inp in enumerate((n, m)):
        eq(apply(ls, rhos[j], 2), (p0, p1)[j], 'perfect state discrimination')
        eq(S.trace(e * rhos[j]), 1 - j, 'effect expectation')
        for a in range(inp):
            for b in range(inp):
                x = unit(inp, inp, a, b)
                y = apply(families[j], x, c)
                eq(apply(ls, y, 2), S.trace(x) * (p0, p1)[j], 'discriminator factorization on full algebra')
    u = S.eye(c) - 2 * qs[0]
    eq(u.H * u, S.eye(c), 'phase unitary Kraus completeness')
    for inp, ks, q in zip((n, m), families, qs):
        for a in range(inp):
            for b in range(inp):
                y = apply(ks, unit(inp, inp, a, b), c)
                eq(q * y * q, y, 'support for every complex matrix unit')
                eq(u * y * u.H, y, 'phase fixes injection on full algebra')
    a, b = v[:, 0], v[:, split]
    cross = a * b.H
    eq(cross * b, a, 'cross-support coherence nonzero')
    require(cross != zero(c), 'nonzero cross-support coherence')
    eq(u * cross * u.H, -cross, 'phase reverses cross coherence')
    require(u * cross * u.H != cross, 'identity and phase conjugation are distinct')
    return c


def zero_cases():
    for n in range(5):
        for z in range(5):
            # A complex-linear map from the 0-dimensional matrix algebra has
            # no free matrix entries; its unique value is zero, with trace 0.
            eq(n * n * 0, 0, 'zero-domain map has zero parameters')
            eq(complete([], 0), S.eye(0), 'zero-source empty Kraus family complete')
            eq(apply([], S.eye(0), z), zero(z), 'unique zero-source channel')
            eq(S.trace(apply([], S.eye(0), z)), S.trace(S.eye(0)), 'zero-source trace preservation')
            for a in range(n):
                for b in range(n):
                    x = unit(n, n, a, b)
                    eq(S.eye(n) * x * S.eye(n), x, 'nonzero leg is identity in zero-input coproduct')
            if n:
                require(complete([], n) != S.eye(n), 'positive source has no zero-target channel')
    print('ZERO CASES: dimensions 0..4, empty Kraus channel, identity leg, positive-to-zero obstruction')


def main():
    fidelity()
    cases = [(n, m, extra, pure) for n in range(1, 4) for m in range(1, 4)
             for extra in (0, 1) for pure in (False, True)]
    cases += [(n, m, extra, pure) for n, m in [(3, 4), (5, 12)]
              for extra in (0, 1) for pure in (False, True)]
    cases += [(3, 4, 3, True), (5, 12, 11, True)]
    dims = [obstruction(*case) for case in cases]
    print(f'PHASE/SUPPORT: {len(cases)} exact channel pairs; dimensions include (1,1), (3,4), (5,12)')
    print(f'PHASE/SUPPORT: candidate vertex sizes {sorted(set(dims))}; full-algebra matrix-unit checks')
    print('PHASE/SUPPORT: rational complex rotated supports, mixed-unitary and pure preparation families')
    zero_cases()
    print(f'EXACT_ASSERTIONS: {checks}')
    print('SCOPE: finite families and proof identities; no enumeration of all CPTP injections')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
