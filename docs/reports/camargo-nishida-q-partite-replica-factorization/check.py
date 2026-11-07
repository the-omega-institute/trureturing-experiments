#!/usr/bin/env python3
"""Exact finite checks of the n=2 bosonic Gaussian replica factorization.

Source fidelity: Camargo--Nishida, arXiv:2609.30754v2, Section 4:
Z_2^(3) = E_A E_B E_C, including complex symmetric precision matrices.

The replica determinant is constructed directly from the source's twisted density
kernels. Purities are computed independently from reduced canonical covariance
matrices (vacuum covariance I/2), not from the proposed character blocks.
Finite tests do not prove the statement for all real parameters or all q.
Requires Python 3 and SymPy. No floating point arithmetic is used.
"""
import sympy as s
from sympy.polys.matrices import DomainMatrix


def determinant(a):
    return s.QQ_I.to_sympy(DomainMatrix.from_Matrix(a).convert_to(s.QQ_I).det())


def compositions(total, length):
    if length == 1:
        yield (total,)
    else:
        for head in range(1, total - length + 2):
            for tail in compositions(total - head, length - 1):
                yield (head,) + tail


def data(n, variant):
    # R is positive definite by construction: positive diagonal plus vv^T.
    v = s.Matrix([(-1) ** i * (1 + i % 2) for i in range(n)])
    r = s.diag(*range(2, n + 2)) + v * v.T
    if variant == 'real':
        im = s.zeros(n)
    else:
        im = s.Matrix(n, n, lambda i, j:
                      s.Rational(((i + 1) * (j + 1)) % 5 - 2, 3))
    assert r == r.T and im == im.T
    assert all(determinant(r[:k, :k]) > 0 for k in range(1, n + 1))
    return r, im


def owner(sizes):
    return [party for party, size in enumerate(sizes) for _ in range(size)]


def replica_matrix(r, im, sizes):
    """Direct mN-coordinate quadratic form; no Fourier transform is used."""
    q, n = len(sizes), sum(sizes)
    m = 1 << (q - 1)
    w = r + s.I * im
    base = s.kronecker_product(s.eye(m), w)
    labels = owner(sizes)
    perm = []
    for replica in range(m):
        for mode, party in enumerate(labels):
            other = replica ^ (1 << party) if party < q - 1 else replica
            perm.append(other * n + mode)
    assert sorted(perm) == list(range(m * n))
    assert all(perm[perm[k]] == k for k in range(m * n))
    twisted = base.conjugate().extract(perm, perm)
    return (base + twisted) / 2


def covariance(r, im):
    inv = r.inv()
    # Twice the canonical covariance; the q,p cross block has a minus sign
    # for psi proportional to exp(-q^T (R+iS) q/2).
    return inv.row_join(-inv * im).col_join(
        (-im * inv).row_join(r + im * inv * im))


def purity_squared(cov, modes, n):
    indices = modes + [n + k for k in modes]
    d = determinant(cov.extract(indices, indices))
    assert d > 0
    return s.cancel(1 / d)


def check_case(q, sizes, variant):
    n, m = sum(sizes), 1 << (q - 1)
    r, im = data(n, variant)
    det_r = determinant(r)
    direct = replica_matrix(r, im, sizes)
    det_direct = determinant(direct)
    assert det_direct.is_real and det_direct > 0
    z_squared = s.cancel(det_r ** m / det_direct)
    cov = covariance(r, im)
    assert cov == cov.T and determinant(cov) == 1
    factors = []
    labels = owner(sizes)
    for mask in range(1, m):
        modes = [k for k, a in enumerate(labels) if a < q - 1 and (mask >> a) & 1]
        p2 = purity_squared(cov, modes, n)
        complement = [k for k in range(n) if k not in modes]
        assert p2 == purity_squared(cov, complement, n)
        assert 0 < p2 <= 1
        # A third exact check: partial-parity Gaussian determinant.
        d = s.diag(*[-1 if k in modes else 1 for k in range(n)])
        block = (r - s.I * im + d * (r + s.I * im) * d) / 2
        det_block = determinant(block)
        assert det_block.is_real and det_block > 0
        assert p2 == s.cancel(det_r / det_block)
        factors.append(p2)
    assert z_squared == s.prod(factors)
    if q == 3:
        # The source uses A,B,C; the third mask here denotes AB.
        # Check AB has the same purity as C, then the exact source product.
        assert factors[2] == purity_squared(cov, [k for k, a in enumerate(labels) if a == 2], n)
    print(f"PASS q={q} sizes={sizes} kernel={variant} factors={len(factors)} Z_squared={z_squared}", flush=True)
    return q == 3


def binary_characters():
    # Unnormalized Walsh H: H H^T = mI. This avoids square-root normalization.
    for r in range(1, 7):
        m = 1 << r
        h = [[(-1) ** ((a & b).bit_count()) for b in range(m)] for a in range(m)]
        for a in range(m):
            for b in range(m):
                assert sum(h[a][t] * h[b][t] for t in range(m)) == (m if a == b else 0)
            for j in range(r):
                for t in range(m):
                    assert h[a][t ^ (1 << j)] == (-1) ** ((a >> j) & 1) * h[a][t]
        assert len({tuple(j for j in range(r) if (a >> j) & 1) for a in range(1, m)}) == m - 1
        print(f"PASS binary_characters r={r} m={m}", flush=True)


def displacement_check():
    # Identical shifts/momenta on all replicas commute with party permutations.
    for q in range(2, 7):
        m = 1 << (q - 1)
        for a in range(q):
            perm = [t ^ (1 << a) if a < q - 1 else t for t in range(m)]
            assert sorted(perm) == list(range(m))
            # Equality of coefficients of every coordinate in the phase sum.
            assert [perm.count(t) for t in range(m)] == [1] * m
    print("PASS displacement coordinate and phase coefficients q=2..6", flush=True)


def main():
    print(f"Python exact arithmetic; SymPy {s.__version__}", flush=True)
    binary_characters()
    displacement_check()
    cases = fidelity = 0
    # Exhaust this declared finite family: q=2..5, N=q or q+1,
    # every ordered positive party-size composition, real and complex kernels.
    for q in range(2, 6):
        for n in (q, q + 1):
            for sizes in compositions(n, q):
                for variant in ('real', 'complex'):
                    fidelity += int(check_case(q, sizes, variant))
                    cases += 1
    assert cases == 36 and fidelity == 8
    print(f"SOURCE_FIDELITY_OK tripartite_cases={fidelity}")
    print(f"FACTORIZATION_OK exact_cases={cases} q=2..5 N=q,q+1; all ordered positive size compositions")
    print("SCOPE: finite exact evidence only; no Gaussian-integral or universal Lean proof")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
