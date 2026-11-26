"""Interval inertia of both SU(5) Weyl kernels on a finite five-periodic torus."""
from collections import Counter
from itertools import combinations, product
from flint import acb, acb_mat, arb, arb_poly, ctx

ROWS = [(0, 2, 1, 4, 3), (3, 4, 2, 1, 0), (4, 3, 1, 2, 0), (2, 4, 1, 0, 3)]
POWERS = [(1, 3, 2, 0, 2), (1, 1, 0, 2, 2), (0, 0, 0, 1, 1), (0, 3, 0, 1, 0)]


def adjoint(a):
    return a.conjugate().transpose()


def identity(n):
    return acb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def tensor(a, b):
    return acb_mat([[a[i, j]*b[k, l] for j in range(a.ncols()) for l in range(b.ncols())]
                    for i in range(a.nrows()) for k in range(b.nrows())])


def field(rows=ROWS, exponents=POWERS):
    values = []
    for permutation, powers in zip(rows, exponents):
        u = acb_mat(5, 5)
        for i, (j, phase) in enumerate(zip(permutation, powers)):
            u[i, j] = acb(0, 1)**phase
        assert u.det() == 1 and u*adjoint(u) == identity(5)
        values.append(u)
    return values+[adjoint(values[0]*values[1]*values[2]*values[3])]


def wedge(u):
    pairs = list(combinations(range(5), 2))
    return acb_mat([[u[i, k]*u[j, l]-u[i, l]*u[j, k] for k, l in pairs] for i, j in pairs])


def stencil(values):
    one = identity(2)
    x, y = acb_mat([[0, 1], [1, 0]]), acb_mat([[0, acb(0, -1)], [acb(0, 1), 0]])
    z = acb_mat([[1, 0], [0, -1]])
    gamma = [tensor(x, one), tensor(y, one), tensor(z, x), tensor(z, y)]
    chirality = gamma[0]*gamma[1]*gamma[2]*gamma[3]
    total = sum(gamma[1:], gamma[0])
    offset = (1-1/arb(5).sqrt())/4
    vectors = [g-total*offset for g in gamma]
    vectors.append(-sum(vectors[1:], vectors[0]))
    # The conjugate transpose of the forward coefficient gives the backward one.
    zero = tensor(chirality, identity(values[0].nrows()))*(arb(61)/40)
    forward = [tensor(chirality*(v-identity(4))/2, u) for v, u in zip(vectors, values)]
    return zero, forward


def positive_roots(poly, cut):
    # With only real nonzero roots, Descartes' sign rule gives the exact count.
    shifted = poly(arb_poly([cut, 1]))
    coeffs = shifted.coeffs()
    assert all(not x.contains(0) for x in coeffs)
    return sum((a > 0) != (b > 0) for a, b in zip(coeffs, coeffs[1:]))


def fibers(values):
    zero, forward = stencil(values)
    c = (arb(5).sqrt()-1)/4
    zeta = acb(c, (1-c*c).sqrt())
    phased = [[matrix*zeta**k for k in range(5)] for matrix in forward]
    for p in product(range(5), repeat=4):
        power = (*p, -sum(p) % 5)
        f = sum((phased[i][power[i]] for i in range(1, 5)), phased[0][p[0]])
        yield zero+f+adjoint(f)


def fiber_indices(values):
    stats = Counter()
    cut = arb(1)/2000
    for hermitian in fibers(values):
        dim = hermitian.nrows()
        poly = arb_poly([x.real for x in hermitian.charpoly().coeffs()])
        above = positive_roots(poly, cut)
        assert above == positive_roots(poly, -cut)
        stats[dim//2-above] += 1
    return stats


def verify():
    ctx.prec = 320
    values = field()
    conjugate = fiber_indices([u.conjugate() for u in values])
    antisymmetric = fiber_indices([wedge(u) for u in values])
    assert conjugate == {-1: 1, 0: 621, 1: 3}
    assert antisymmetric == {-1: 6, 0: 606, 1: 13}
    indices = tuple(sum(j*n for j, n in counts.items()) for counts in (conjugate, antisymmetric))
    assert indices == (2, 7) and (-indices[0]+2*indices[1]) % 5 == 2
    print("Certified fiber counts, index pair and determinant-line center charge:",
          conjugate, antisymmetric, indices, 2)


if __name__ == "__main__":
    verify()
