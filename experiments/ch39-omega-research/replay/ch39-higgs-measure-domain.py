"""Exact examples for the link correction and finite Higgs-fiber obstruction."""
import sympy as sp


def phase(t):
    return sp.cancel((1 + sp.I*t) / (1 - sp.I*t))


def correction(v, w):
    identity = sp.eye(5)
    inner = (v.H*w)[0].expand(complex=True).simplify()
    radius = sp.Abs(inner)
    unit = sp.simplify(inner / radius)
    # The positive-real-part square root of a unit phase in the right half-plane.
    def unit_sqrt(z):
        return sp.simplify((1 + z) / sp.sqrt(2 + z + sp.conjugate(z)))

    fourth = unit_sqrt(unit_sqrt(sp.conjugate(unit)))
    direction = v*v.H
    B = (unit*direction + fourth*(identity - direction)).applyfunc(sp.simplify)
    u = unit*v
    K = w*u.H - u*w.H
    R = (identity + K + K*K/(1 + radius)).applyfunc(sp.simplify)
    return (R*B).applyfunc(sp.simplify)


def link_examples():
    basis = sp.eye(5)
    v = basis[:, 0]
    root = phase(sp.Rational(1, 32))
    twist = sp.expand_complex(root**4).simplify()
    c, s = sp.Rational(24, 25), sp.Rational(7, 25)
    frame = sp.eye(5)
    frame[:2, :2] = sp.Matrix([[c, sp.I*s], [sp.I*s, c]])
    for unit, cosine, sine in ((1, 1, 0), (twist, 1, 0), (1, c, s), (twist, c, s)):
        w = unit*cosine*v + sp.I*sine*basis[:, 2]
        S = correction(v, w)
        assert (S*v - w).applyfunc(sp.simplify) == sp.zeros(5, 1)
        assert (S.H*S).applyfunc(sp.simplify) == basis and sp.simplify(S.det()) == 1
        # Exact operator-norm bound via positive semidefiniteness.
        eta2 = sp.simplify(((v-w).H*(v-w))[0])
        difference = S - basis
        assert (3*eta2*basis - difference.H*difference).applyfunc(sp.simplify).is_positive_semidefinite
        if unit == twist and sine == s:
            transformed = correction(frame*v, frame*w)
            assert (transformed - frame*S*frame.H).applyfunc(sp.simplify) == sp.zeros(5)


def projection_and_fiber():
    c, s = sp.Rational(24, 25), sp.Rational(7, 25)
    v = sp.Matrix([1, 0])
    w = sp.Matrix([c, s])
    P, Q = v*v.H, w*w.H
    delta = 2*(P-Q)
    assert (sp.eye(2) - delta**2).is_positive_definite
    assert Q*v/sp.sqrt((v.H*Q*v)[0]) == w
    T = sp.symbols("T")
    for order in (2, 5):
        unity = sp.cyclotomic_poly(order, T)
        for charge in range(1, order):
            assert sp.gcd(unity, T**charge - 1) == 1
    # A two-column edge residual with nonorthogonal columns.
    residual = sp.Matrix([[c-1, sp.I*s], [sp.I*s, c-1], [s, s]])
    gram = residual.H*residual
    assert max(gram.eigenvals()) == (c-1)**2 + 3*s*s


if __name__ == "__main__":
    link_examples()
    projection_and_fiber()
    print("SU(5) link bounds, equivariance, spectral transport and cyclic charge checks passed.")
