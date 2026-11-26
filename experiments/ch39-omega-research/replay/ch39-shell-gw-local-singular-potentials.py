"""Exclude flavor mixing for independent homogeneous local potentials of singular values."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb_mat, arb

folder = Path(__file__).resolve().parent
vac = run_path(str(folder/"ch39-shell-gw-purity-vacuum.py"))
four = run_path(str(folder/"ch39-shell-gw-four-mode-obstruction.py"))
m = vac["m"]
k, g = vac["k"], m["g"]
adj, hs, rows, count = (m[x] for x in ("adj","hs","rows","count"))
as_array = four["float_matrix"]


def transported_cartan_cost():
    _, phases = vac["ground"]()
    phase = np.repeat([complex(z) for z in phases],8)
    kf = as_array(k)*phase.conj()[:,None]*phase[None,:]
    _, vectors = np.linalg.eigh(kf.real)
    generators = np.array([as_array(v) for v in g])

    def trial(angle):
        coeff = vectors[:,0]*np.cos(angle)+vectors[:,1]*np.sin(angle)
        matrices = np.einsum("va,aij->vij", coeff.reshape(6,8),generators)
        sites = np.linalg.eigh(matrices)[1]
        fields = sites[:,None] @ generators[:2] @ sites[:,None].conj().transpose(0,1,3,2)
        coefficients = fields.reshape(6,2,9)@generators.reshape(8,9).conj().T
        q = coefficients.transpose(1,0,2).reshape(2,48).T
        return np.linalg.eigvalsh((q.conj().T@kf@q).real/6)[-1], sites

    _, sites = min((trial(theta) for theta in np.linspace(0,2*np.pi,91)),key=lambda q:q[0])
    unitaries = [four["exact_frame"](site) for site in sites]
    cols = []
    for d in g[:2]:
        col = []
        for phase, u in zip(phases,unitaries):
            block = u*d*adj(u)
            col.extend([phase*hs(b,block)] for b in g)
        cols.append(acb_mat(col))
    gram = acb_mat([[hs(x,k*y).real/6 for y in cols] for x in cols])
    assert count(rows(gram), arb(".015908")) == 2


def low_gain():
    sources = [(v,0) for v in m["frame"](m["sing"],arb(".02"))]
    sources += [(v,1) for v in m["frame"](k,arb(".02"))]
    left,right = (m["yuk"]["zero_basis"](s) for s in (1,-1))
    bases = [[m["yuk"]["identity"](3)/arb(3).sqrt()], g]
    phases = m["phases"]
    matrices = []
    for v,kind in sources:
        bs = bases[kind]
        y = acb_mat(3,3)
        for i,phase in enumerate(phases):
            h = sum((v[len(bs)*i+j,0]*b for j,b in enumerate(bs)),acb_mat(3,3))
            l = acb_mat([[left[3*i+t,j] for j in range(3)] for t in range(3)])
            r = acb_mat([[right[3*i+t,j] for j in range(3)] for t in range(3)])
            y += phase*l.transpose()*h.transpose()*r
        matrices.append(y)
    metric = acb_mat([[hs(a,b) if i==j else 0 for b,j in sources] for a,i in sources])
    gram = acb_mat([[hs(a,b) for b in matrices] for a in matrices])
    assert count(rows(gram-arb(".40")**2*metric),arb(0)) == 0


if __name__ == "__main__":
    transported_cartan_cost()
    low_gain()
    assert count(rows(k),arb("1.0832")) == 2
    assert count(rows(m["sing"]),arb("1.0832")) == 1
    # The next modes in either representation cost at least this excess above
    # their respective ground energies.
    gap = arb("1.0832")-arb(".01588064")
    fraction = (arb(2)/3*(arb(".015908")-arb(".01588062"))/gap).sqrt()
    rel = fraction/(arb(".40")*(1-fraction**2).sqrt()-fraction)
    assert fraction < arb(".00414") and rel < arb(".0105")
    print("high scalar fraction <",fraction.upper().str(8),
          "high Yukawa fraction <",rel.upper().str(8),"required > .01437")
