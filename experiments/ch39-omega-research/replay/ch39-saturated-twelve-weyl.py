"""Ball-certified ranks on a 2^4 box and a nonzero locally saturated insertion."""
from itertools import product
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, ctx

BASE = run_path(str(Path(__file__).with_name("ch39-rough-overlap-index.py")))
ADJ = BASE["adjoint"]
POS = list(product(range(2), repeat=4))


def shifts(a,b):
    steps = [acb_mat(16,16) for _ in range(4)]
    for i, x in enumerate(POS):
        exponents = (-2*x[0]*x[1], x[0], -2*x[2]*x[3], x[2])
        for mu in range(4):
            y = list(x)
            y[mu] = 1-y[mu]
            power = (a if mu < 2 else b)*exponents[mu]
            steps[mu][i,POS.index(tuple(y))] = (-1)**x[mu]*acb(0,1)**power
    steps.append(ADJ(steps[0]*steps[1]*steps[2]*steps[3]))
    return steps


def norm_bound(matrix):
    n = matrix.nrows()
    return n*max(abs(matrix[i,j]).upper() for i in range(n) for j in range(n))


def certified_frame(h):
    n = h.nrows()
    approx = np.array([[complex(h[i,j]) for j in range(n)] for i in range(n)])
    vals, vecs = np.linalg.eigh(approx)
    scale = 2**44
    raw = np.rint(vals*scale).astype(np.int64)
    re = np.rint(vecs.real*scale).astype(np.int64)
    im = np.rint(vecs.imag*scale).astype(np.int64)
    v = acb_mat([[acb(int(re[i,j]),int(im[i,j])) for j in range(n)] for i in range(n)])/scale
    diag = acb_mat(n,n)
    for i in range(n):
        diag[i,i] = arb(int(raw[i]))/scale
    one = BASE["identity"](n)
    assert max(abs(int(y)) for y in raw) < 16*scale
    assert norm_bound(ADJ(v)*v-one) < arb(2)**-23
    assert norm_bound(h-v*diag*ADJ(v)) < arb(2)**-23
    return [arb(int(y))/scale for y in raw], v


def verify():
    ctx.prec = 110
    one = BASE["identity"](64)
    chirality = acb_mat(64,64)
    for i in range(64):
        chirality[i,i] = (-1,1,1,-1)[i//16]
    for a,b in ((0,0),(1,0),(0,1),(1,1)):
        zero, forward = BASE["stencil"](shifts(a,b))
        f = sum(forward[1:],forward[0])
        h = zero*100/61+f+ADJ(f)  # 5-m=5/2 when m=5/2.
        eigen, v = certified_frame(h)
        assert min(abs(e) for e in eigen) > arb(".041")
        index = 32-sum(e>0 for e in eigen)
        expected_nullity = int((a,b)==(1,1))
        assert index == -expected_nullity
        signs = acb_mat(64,64)
        for i,e in enumerate(eigen):
            signs[i,i] = 1 if e > 0 else -1
        approx_d = one+chirality*v*signs*ADJ(v)
        squared, _ = certified_frame(ADJ(approx_d)*approx_d)
        assert sum(t <= 1 for t in squared) == expected_nullity
        print((a,b), "I =",index,"positive squared singular values > .99")


if __name__ == "__main__":
    verify()
