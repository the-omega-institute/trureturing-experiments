"""Exact adjoint-spurion realization of the shell trace on every SU(5) gauge direction."""
from itertools import combinations
import sympy as sp

P = sp.diag(0, 0, 0, 1, 1)
I = sp.eye(5)
Z = I-2*P
H = 2*I-5*P
PAIRS = list(combinations(range(5), 2))
SIGNS = [Z[i, i] for i in range(5)]
Z_W = sp.diag(*SIGNS, *(SIGNS[i]*SIGNS[j] for i, j in PAIRS))
T = sp.Rational(1, 15)


def weyl(x):
    def entry(i, j):
        a, b = PAIRS[i]
        c, d = PAIRS[j]
        return x[a, c]*int(b == d)-x[b, c]*int(a == d)+int(a == c)*x[b, d]-int(b == c)*x[a, d]
    return sp.diag(-x.T, sp.Matrix(10, 10, entry))


def tr_product(x, y):
    return sum(x[i, j]*y[j, i] for i, j in x.todok())


def verify():
    basis = []
    sectors = []
    for i, j in PAIRS:
        for value in (1, sp.I):
            x = sp.zeros(5)
            x[i, j], x[j, i] = value, sp.conjugate(value)
            basis.append(x)
            sectors.append('color' if j < 3 else 'weak' if i >= 3 else 'broken')
    for i in (0, 1, 3):
        x = sp.zeros(5)
        x[i, i], x[i+1, i+1] = 1, -1
        basis.append(x)
        sectors.append('weak' if i == 3 else 'color')
    basis.append(H)
    sectors.append('hypercharge')
    assert len(basis) == 24 and sp.Matrix.hstack(*(x.reshape(25, 1) for x in basis)).rank() == 24
    assert Z.det() == 1 and Z_W.trace() == -1
    assert (H-2*I)*(H+3*I) == sp.zeros(5)
    charged = [weyl(x) for x in basis]
    blocks = [P*x*P for x in basis]
    shift = {'hypercharge': -sp.Rational(3, 5), 'weak': 1, 'color': 0, 'broken': 0}
    ratios = {sector: 1+T*d for sector, d in shift.items()}
    for i, x in enumerate(basis):
        for j in range(i, len(basis)):
            base = tr_product(x, basis[j])
            if sectors[i] != sectors[j]:
                assert base == 0
            native = tr_product(charged[i], charged[j])
            # Symmetric bilinear trace, including broken directions which need not commute with Z_W.
            weighted = (tr_product(Z_W*charged[i], charged[j])+tr_product(Z_W*charged[j], charged[i]))/2
            correction = tr_product(blocks[i], blocks[j])-blocks[i].trace()*blocks[j].trace()
            assert native == 4*base and weighted == -4*correction
            assert correction == shift[sectors[i]]*base
    assert ratios == {'hypercharge': sp.Rational(24, 25), 'weak': sp.Rational(16, 15), 'color': 1, 'broken': 1}
    weak_angle = 3*ratios['weak']/(5*ratios['hypercharge']+3*ratios['weak'])
    assert weak_angle == sp.Rational(2, 5)
    # H insertion multipliers of color, weak and broken generators under the fundamental trace.
    multipliers = {}
    for sector in ('color', 'weak', 'broken'):
        x = basis[sectors.index(sector)]
        multipliers[sector] = (H*x*x).trace()/(x*x).trace()
    assert multipliers == {'color': 2, 'weak': -3, 'broken': -sp.Rational(1, 2)}
    print('300 entries of the symmetric SU(5) kinetic form; eigenvalue ratios', ratios, 'sin^2(theta_W)', weak_angle)


if __name__ == '__main__':
    verify()
