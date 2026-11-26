"""Order-five electromagnetic flux on a 5^4 torus and the weak Cartan determinant charge."""
from itertools import combinations, combinations_with_replacement, product
from pathlib import Path
from runpy import run_path
from flint import arb, acb, acb_mat, arb_poly, ctx

base = run_path(str(Path(__file__).with_name("ch39-rough-overlap-index.py")))
whole = run_path(str(Path(__file__).with_name("ch39-uniform-center-obstruction.py")))
INDICES = (0, 2, 1, 1, 2)


def magnetic_irrep(charge):
    eye = base["identity"](5)
    c = (arb(5).sqrt()-1)/4
    z = acb(c, (1-c*c).sqrt())
    shift, clock = acb_mat(5, 5), acb_mat(5, 5)
    for i in range(5):
        shift[i, (i+1) % 5] = 1
        clock[i, i] = z**((charge*i) % 5)
    tensor = base["tensor"]
    u = [tensor(shift, eye), tensor(clock, eye), tensor(eye, shift), tensor(eye, clock)]
    return u+[base["adjoint"](u[0]*u[1]*u[2]*u[3])]


def charge_index(charge):
    u = magnetic_irrep(charge)
    mass, forward = base["stencil"](u)
    f = sum(forward[1:], forward[0])
    h = mass+f+base["adjoint"](f)
    chirality = mass*(arb(40)/61)
    det = arb_poly([v.real for v in (chirality*h).charpoly().coeffs()])
    whole["exclude_zero"](det, -arb(21)/1000, arb(25)/1000)
    poly = arb_poly([v.real for v in h.charpoly().coeffs()])
    cut = arb(1)/100
    npos = base["positive_roots"](poly, cut)
    assert npos == base["positive_roots"](poly, -cut)
    return 50-npos


def stabilizer_weights(charges):
    return tuple(-INDICES[-x % 5]+sum(INDICES[(x+y) % 5] for j, y in enumerate(charges) if i != j)
                 for i, x in enumerate(charges))


def verify():
    ctx.prec = 896
    for charge in range(1, 5):
        index = charge_index(charge)
        assert index == INDICES[charge]
        print("Charge and certified irreducible magnetic index:", charge, index, flush=True)
    classes = [p for p in combinations_with_replacement(range(5), 5) if sum(p) % 5 == 0]
    assert len(classes) == 26 and sum(len(set(stabilizer_weights(p))) == 1 for p in classes) == 11
    ordered = [(*p, -sum(p) % 5) for p in product(range(5), repeat=4)]
    assert sum(len(set(stabilizer_weights(p))) > 1 for p in ordered) == 460
    photon = (1, 1, 1, -3, 0)
    v = stabilizer_weights(photon)
    assert v == (3, 3, 3, 3, 7)
    i5 = 25*sum(INDICES[-q % 5] for q in photon)
    i10 = 25*sum(INDICES[(p+q) % 5] for p, q in combinations(photon, 2))
    assert (i5, i10) == (175, 325) and (-i5+2*i10) % 5 == 0
    for generator, total in [((0, 0, 0, 1, -1), -100), ((-2, -2, -2, 3, 3), 300), (photon, 0)]:
        assert 25*sum(a*b for a, b in zip(v, generator)) == total
    assert arb(3454)/1000 < (25-5*arb(5).sqrt())/4


if __name__ == "__main__":
    verify()
