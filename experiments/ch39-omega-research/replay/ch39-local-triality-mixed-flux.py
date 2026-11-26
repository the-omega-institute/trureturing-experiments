"""Certified mixed continuous and local shell-triality flux on a five-periodic box."""
from itertools import product
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, arb_poly, ctx

base = run_path(str(Path(__file__).with_name("ch39-rough-overlap-index.py")))
FIRST = (0, 0, 1, 1, 2)
SECOND = (0, 1, 1, 2, 2)
EXPECTED = (0,1,1,1,1, 1,1,1,1,2, 0,0,0,0,0, 1,0,1,0,0, 1,1,0,1,1)


def phase_fibers(first=FIRST, second=SECOND):
    five = base["identity"](5)
    shift = acb_mat(5,5)
    for i in range(5):
        shift[i,(i+1)%5] = 1
    omega = acb(-arb(1)/2, arb(3).sqrt()/2)
    clocks = []
    for pattern in (first, second):
        clock = acb_mat(5,5)
        for i, exponent in enumerate(pattern):
            clock[i,i] = omega**exponent
        clocks.append(clock)
    tensor = base["tensor"]
    steps = [tensor(shift,five), tensor(clocks[0],five),
             tensor(five,shift), tensor(five,clocks[1])]
    steps.append(base["adjoint"](steps[0]*steps[1]*steps[2]*steps[3]))
    mass, forward = base["stencil"](steps)
    real = (arb(5).sqrt()-1)/4
    root = acb(real, (1-real**2).sqrt())
    phased = {i:[forward[i]*root**k for k in range(5)] for i in (1,3,4)}
    for j,k in product(range(5),repeat=2):
        f = forward[0]+phased[1][j]+forward[2]+phased[3][k]+phased[4][(-j-k)%5]
        yield mass+f+base["adjoint"](f)


def verify():
    ctx.prec = 896
    cut = arb(1)/10**6
    indices = []
    for matrix in phase_fibers():
        poly = arb_poly([c.real for c in matrix.charpoly().coeffs()])
        above = base["positive_roots"](poly, cut)
        assert above == base["positive_roots"](poly, -cut)
        indices.append(50-above)
    assert tuple(indices) == EXPECTED and sum(indices) == 16
    jacobian = 2*sum(indices)
    assert jacobian % 3 == 2
    print("Certified 25 magnetic fibers, gap > 1e-6, I(S)=16 and local triality exponent 32.")


if __name__ == "__main__":
    verify()
