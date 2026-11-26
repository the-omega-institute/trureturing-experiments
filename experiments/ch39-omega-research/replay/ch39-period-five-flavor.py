"""Period-five Yukawa maps on the central momentum orbit of the simplex overlap kernel."""
import sympy as sp
from sympy.polys.matrices import DomainMatrix

FIELD = sp.QQ.cyclotomic_field(5)
ROOT = FIELD.unit
LIGHT = {0, 1, 4}


def rank_on(sources, powers):
    rows = [[ROOT**((s*j) % 5) for j in powers] for s in sorted(sources)]
    return DomainMatrix(rows, (len(rows), len(powers)), FIELD).rank()


def study_stencil(powers):
    compressed = []
    isolated = []
    for transfer in range(5):
        active = {s for s in LIGHT if (s+transfer) % 5 in LIGHT}
        cross = {s for s in range(5) if (s in LIGHT) != ((s+transfer) % 5 in LIGHT)}
        compressed.append(rank_on(active, powers))
        isolated.append(rank_on(active | cross, powers)-rank_on(cross, powers))
    return compressed, isolated


def rank_certificates():
    expected = [
        ([0], [1, 1, 1, 1, 1], [1, 0, 0, 0, 0]),
        ([0, 1], [2, 2, 1, 1, 2], [2, 0, 0, 0, 0]),
        ([-1, 0, 1], [3, 2, 1, 1, 2], [3, 1, 0, 0, 1]),
        ([-2, -1, 0, 1, 2], [3, 2, 1, 1, 2], [3, 2, 1, 1, 2]),
    ]
    for powers, free, isolated in expected:
        assert study_stencil(powers) == (free, isolated)
        print("Powers, compressed dimension, isolated dimension:", powers, sum(free), sum(isolated))
    # Exact inverse transform yields all five source projectors.
    for j in range(5):
        for k in range(5):
            numerator = sum((ROOT**((h*(j-k)) % 5) for h in range(-2, 3)), FIELD.zero)
            assert numerator == FIELD.convert(5*int(j == k))


def light_matrix(coeff):
    # Both legs are normalized by the inverse square root of the free kinetic coefficient.
    a = sp.Symbol("a", positive=True)
    normalizer = sp.diag(1, a, 1)
    fourier = sp.Matrix(3, 3, lambda i, j: coeff[(i-j) % 5])
    return normalizer*fourier*normalizer


def flavor_certificates():
    a = sp.Symbol("a", positive=True)
    up = light_matrix([1, 1, 0, 0, 0])
    down = light_matrix([1, 1, sp.I, 0, 0])
    h_up, h_down = up*up.H, down*down.H
    assert sp.factor((h_up*h_down-h_down*h_up).det()) == -2*sp.I*a**10
    assert up.det() == down.det() == a**2
    coeff = sp.symbols("c:5")
    mirror = sp.eye(3)[::-1, :]
    first = light_matrix(coeff)
    conjugate_profile = light_matrix([sp.conjugate(coeff[-k % 5]) for k in range(5)])
    assert conjugate_profile == first.H == mirror*first.conjugate()*mirror


if __name__ == "__main__":
    rank_certificates()
    flavor_certificates()
    print("The single-profile obstruction and the nonzero CP invariant hold at every positive residue ratio.")
