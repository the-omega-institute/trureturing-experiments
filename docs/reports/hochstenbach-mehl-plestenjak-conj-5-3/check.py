#!/usr/bin/env python3
"""Exact finite checks for arXiv:2406.07109v1, Conjecture 5.3.

Run with Python 3 and SymPy. All entries and calculations are exact.
The finite grid is evidence for a proposed identity, not a universal proof.
"""
from itertools import product
import hashlib
from pathlib import Path
import sys
import sympy as sp

x, t = sp.symbols('lambda tau')


def det(a):
    return sp.expand(a.det(method='domain-ge'))


def equal(a, b):
    assert sp.Poly(sp.expand(a - b), x, t).is_zero, (a, b)


def source_fidelity():
    # Section 5, displayed 3-by-3 pencil: right minimal index exactly one.
    p = sp.Matrix([[0, 1, -x], [1, 0, 0], [-x, 0, 0]])
    a, b = p.subs(x, 0), -p.diff(x)
    kernel = sp.Matrix([0, x, 1])
    assert p * kernel == sp.zeros(3, 1)
    assert det(p[:2, :2]) == -1  # normal rank is two
    assert a.col_join(b).rank() == 3  # no nonzero constant kernel vector
    u = sp.ones(3, 1)
    q = p + (2 - x) * u * u.T
    equal(det(q), (x - 2) * (x + 1)**2)
    assert q.subs(x, -1).rank() == 2
    assert (u.T * q.subs(x, -1).nullspace()[0])[0] == 0
    assert det(q).subs(x, 2) == 0
    print('FIDELITY: source 3x3 pencil has normal rank 2, minimal index 1;')
    print('  prescribed eigenvalue 2; random eigenvalue -1 is double, geometric multiplicity 1.')
    print('SOURCE_TEX_SIGN_DISCREPANCY: the printed example says +1, but the displayed pencil gives -1.')
    # The next real rank-one example has the same sign discrepancy.
    p5 = sp.Matrix([[0, 0, 1, -x, 0], [0, 0, 0, 1, -x],
                    [1, 0, 0, 0, 0], [-x, 1, 0, 0, 0], [0, -x, 0, 0, 0]])
    equal(det(p5 + (2 - x) * sp.ones(5, 1) * sp.ones(1, 5)),
          (2 - x) * (x*x + x + 1)**2)
    print('SOURCE_TEX_SIGN_DISCREPANCY: the printed 5x5 example says real part +1/2; exact factor is lambda^2+lambda+1.')


def canonical(indices, r):
    """R followed by short then long coordinates in each symmetric singular block."""
    blocks = [r]
    kernel_cols = []
    cursor = r.rows
    n = r.rows + sum(2*m + 1 for m in indices)
    for m in indices:
        ell = sp.zeros(m, m + 1)
        for i in range(m):
            ell[i, i] = -x
            ell[i, i+1] = 1
        block = sp.zeros(2*m+1)
        block[:m, m:] = ell
        block[m:, :m] = ell.T
        blocks.append(block)
        col = sp.zeros(n, 1)
        for j in range(m + 1):
            col[cursor + m + j] = x**j
        kernel_cols.append(col)
        cursor += 2*m+1
    return sp.diag(*blocks), sp.Matrix.hstack(*kernel_cols)


def identity(indices, r, d, u):
    p, nker = canonical(indices, r)
    g = det(nker.T * u)
    expected = (-1)**sum(indices) * t**len(indices) * det(r) * det(d) * g**2
    actual = det(p + t*u*d*u.T)
    equal(actual, expected)
    return actual, g


def finite_grid():
    # All ordered lists of length 1..3 with entries 0..2, all r=0..2,
    # subject to total pencil size <=7. Three explicitly generated U values.
    shapes = cases = 0
    for k in range(1, 4):
        for indices in product(range(3), repeat=k):
            for rsize in range(3):
                n = rsize + 2*sum(indices) + k
                if n > 7:
                    continue
                shapes += 1
                r = sp.diag(*[x - (j+2) for j in range(rsize)])
                d = sp.diag(*[x - (j+7) for j in range(k)])
                # Non-diagonal symmetric prescribed pencil, related by an invertible constant congruence.
                v = sp.eye(k)
                for j in range(k-1):
                    v[j, j+1] = j+1
                d = v*d*v.T
                for mode in range(3):
                    if mode == 0:
                        u = sp.zeros(n, k)
                    else:
                        u = sp.Matrix(n, k, lambda i, j: ((i+1)*(j+2) + mode*i*i + j*j) % 7 - 3)
                    identity(indices, r, d, u)
                    cases += 1
    print(f'IDENTITY_GRID: {shapes} block-size choices; {cases} exact polynomial identities in lambda and tau.')
    # A full symbolic U (all 8 independent entries), k=2 and indices (0,1).
    indices = (0, 1)
    p, nker = canonical(indices, sp.zeros(0))
    u = sp.Matrix(4, 2, sp.symbols('u0:8'))
    d = sp.Matrix([[2-x, 1], [1, -x]])
    g = sp.expand((nker.T*u).det())
    equal(det(p+t*u*d*u.T), -t**2 * det(d) * g**2)
    print('SYMBOLIC_U: indices (0,1), non-diagonal D, all 8 U entries indeterminate: residual 0.')


def generic_witness(indices, r):
    p, nker = canonical(indices, r)
    u = sp.zeros(p.rows, len(indices))
    cursor, root = r.rows, 4
    for j, m in enumerate(indices):
        q = sp.prod(x-a for a in range(root, root+m))
        root += m
        for e in range(m+1):
            u[cursor+m+e, j] = sp.Poly(q, x).nth(e)
        cursor += 2*m+1
    # Freely add entries in regular and short coordinates, which do not enter N^T U.
    for j in range(len(indices)):
        for i in range(r.rows):
            u[i, j] = i+j+1
    cursor = r.rows
    for m in indices:
        for i in range(m):
            for j in range(len(indices)):
                u[cursor+i, j] = 2+i-j
        cursor += 2*m+1
    return u


def multiplicity_checks():
    indices = (1, 2)
    r = sp.Matrix([[x-2]])
    d = sp.Matrix([[1-x, 1], [1, -x]])
    u = generic_witness(indices, r)
    actual, g = identity(indices, r, d, u)
    equal(g, (x-4)*(x-5)*(x-6))
    assert sp.discriminant(g, x) == 4
    assert sp.resultant(g, det(r)*det(d), x) != 0
    equal(actual, -t**2*(x-2)*(x*x-x-1)*g**2)
    for tau in (1, 2, -3):
        q = sp.Poly(actual.subs(t, tau), x)
        for root in (4, 5, 6):
            assert q.eval(root) == q.diff().eval(root) == 0
            assert q.diff((x, 2)).eval(root) != 0
    print('GENERIC_WITNESS: n=9 k=2, g=(lambda-4)(lambda-5)(lambda-6), discriminant=4, resultant nonzero.')
    print('  det=-tau^2*(lambda-2)*(lambda^2-lambda-1)*g^2; exact multiplicity 2 at 4,5,6.')
    # R and D each have one eigenvalue at infinity: combined infinity multiplicity two.
    ri = sp.Matrix([[1]])
    di = sp.diag(x-2, 1)
    qi, gi = identity(indices, ri, di, generic_witness(indices, ri))
    assert sp.degree(qi, x) == 7
    assert 9 - sp.degree(qi, x) == (1-sp.degree(det(ri), x))+(2-sp.degree(det(di), x)) == 2
    equal(gi, g)
    print('INFINITY: n=9 determinant degree 7, true infinity multiplicity 1 plus prescribed 1; random roots finite and double.')
    # Repeated and coincident fixed spectra are allowed, and D need not have simple spectrum.
    rr = sp.diag(x-2, x-2)
    dd = sp.diag(x-2, x-2)
    qc, gc = identity((0, 1), rr, dd, generic_witness((0, 1), rr))
    assert sp.Poly(qc, x).degree() == 6
    equal(qc, -t**2*(x-2)**4 * gc**2)
    assert gc.subs(x, 2) != 0
    print('FIXED_COINCIDENCE: two true plus two prescribed copies at lambda=2 give multiplicity 4, additional root remains double.')
    # M=0, n=k: no random eigenvalues, including k=n=1.
    for k in range(1, 4):
        r0 = sp.zeros(0)
        u0 = sp.eye(k)
        q0, g0 = identity((0,)*k, r0, sp.eye(k), u0)
        assert g0 == 1
        equal(q0, t**k)
    print('ZERO_INDICES: M=0 and empty regular part checked for k=n=1,2,3; no random roots.')


def main():
    print('Python', sys.version.split()[0], 'SymPy', sp.__version__)
    print('check.py sha256', hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    source_fidelity()
    finite_grid()
    multiplicity_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
