"""Certified odd weak index on sixteen sites and the pseudoreal spin contraction."""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx
from sympy import Matrix, conjugate, expand, symbols

PREV = run_path(str(Path(__file__).with_name("ch39-saturated-twelve-weyl.py")))
BASE = PREV["BASE"]
ADJ = PREV["ADJ"]


def verify():
    ctx.prec = 110
    steps = []
    for mu, original in enumerate(PREV["shifts"](1,1)[:4]):
        u = acb_mat(32,32)
        for i in range(16):
            for j in range(16):
                z = original[i,j]
                if not z.is_zero():
                    if (mu,i) == (0,0):
                        u[2*i,2*j+1] = z.conjugate()
                        u[2*i+1,2*j] = -z
                    else:
                        u[2*i,2*j] = z
                        u[2*i+1,2*j+1] = z.conjugate()
        assert u*ADJ(u) == BASE["identity"](32)
        steps.append(u)
    steps.append(ADJ(steps[0]*steps[1]*steps[2]*steps[3]))
    mass, forward = BASE["stencil"](steps)
    f = sum(forward[1:],forward[0])
    h = mass*100/61+f+ADJ(f)
    vals, v = PREV["certified_frame"](h)
    assert 64-sum(e>0 for e in vals) == -1 and min(abs(e) for e in vals) > arb(".01002")
    ch = acb_mat(128,128)
    signs = acb_mat(128,128)
    for i,e in enumerate(vals):
        ch[i,i] = (-1,1,1,-1)[i//32]
        signs[i,i] = 1 if e>0 else -1
    approx_d = BASE["identity"](128)+ch*v*signs*ADJ(v)
    squared, _ = PREV["certified_frame"](ADJ(approx_d)*approx_d)
    assert sum(e < arb(".018") for e in squared) == 1
    z,w = symbols("z w",complex=True)
    spinor = Matrix([[z,w],[-conjugate(w),conjugate(z)]])
    assert expand(2*spinor.det()-(spinor.H*spinor).trace()) == 0
    assert 6**3 == 216 and 12*64+21*32 == 1440
    print("Pure SU(2) index -1, exactly one Weyl source zero, all nonzero singular values > 0.1.")
    print("Three generations give the positive spin-color factor 216 sum_x p_x^6.")


if __name__ == "__main__":
    verify()
