"""Finite identities for a gauge-invariant Wilson-Yukawa term."""
from collections import Counter
from itertools import product
from sympy import I, Matrix, Rational as Q, conjugate, diag, eye, simplify, sqrt, zeros


def join(h):
    return Matrix.hstack(Matrix([-conjugate(h[1]),conjugate(h[0])]),h)


def verify():
    w = Matrix([[3,4],[-4,3]])/5 * diag(I,1)
    shell = (-1+sqrt(3)*I)/2
    h = Matrix([1+I, 2-I])/3
    b = join(h)
    assert simplify(b.H*b - (h.H*h)[0]*eye(2)) == zeros(2)
    # The quark columns take sources u_R=C D^-1 S and d_R=C S^-1.
    # det(W)=D^-1. Factor the shared color matrix out.
    d = 1/w.det()
    right = diag(shell/d, 1/shell)
    assert simplify(join(shell*w*h)*right - w*b) == zeros(2)
    # Leptons take sources e_R=D S and nu_R=1.  L is anti-W with S^-1.
    bl = join(h.conjugate())
    assert simplify(join((shell*w*h).conjugate())*diag(d*shell,1)
                    - w.conjugate()*bl/shell) == zeros(2)
    assert b.det() != 0
    corner_masses = Counter(2*sum(bits) for bits in product((0,1), repeat=4))
    assert corner_masses == {0:1, 2:4, 4:6, 6:4, 8:1}
    # On the nonzero-momentum corners the mass is proportional to Higgs length.
    for amplitude in [Q(0), Q(1,3), Q(2)]:
        assert simplify((amplitude*b).H*(amplitude*b)-
                        amplitude**2*(h.H*h)[0]*eye(2)) == zeros(2)
    print("Quark and lepton Higgs maps covariant; corner masses: 0,2,4,6,8 with 1,4,6,4,1.")


if __name__ == "__main__":
    verify()
