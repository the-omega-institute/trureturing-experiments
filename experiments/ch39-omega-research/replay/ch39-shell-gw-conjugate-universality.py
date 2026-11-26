"""Universal simultaneous Yukawa readout from a single covariantly conjugated profile."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

c = run_path(str(Path(__file__).with_name("ch39-shell-gw-conjugate-cut.py")))
m = c["m"]
adj,hs,count,rows = c["adj"],c["hs"],c["count"],c["rows"]
p = c["singlet_phases"]()


def real_matrix(x):
    return acb_mat([[float(v) for v in row] for row in x])


def real_family_kernel(cut,n):
    high = c["frame"](-c["family"],48-n,-cut)
    re, im = [], []
    for col in high:
        z = [p[i//8].conjugate()*col[i,0] for i in range(48)]
        re.append([v.real for v in z])
        im.append([v.imag for v in z])
    M = acb_mat(re+im)
    numeric = c["a"]["four"]["float_matrix"](M).real
    U,s,Vh = np.linalg.svd(numeric, full_matrices=True)
    # Scale the constraint rows before the interval solve.
    good = real_matrix((U/s).T)*M
    gram = good*good.transpose()
    assert count(rows(gram),arb(".99")) == 0
    trial = real_matrix(Vh[len(s):].T)
    return trial-good.transpose()*gram.solve(good*trial)


def singlet_band(cut):
    # Cyclic Fourier and the parity of the two site orbits diagonalize the
    # singlet transport in the ground-state gauge, as proved in the note.
    omega = acb(0,2*arb.pi()/3).exp()
    out = []
    for f in (0,1,-1):
        for sign in (1,-1):
            z = [sign**(i%2)*omega**(f*(i//2)) for i in range(6)]
            v = acb_mat([[p[i]*z[i]] for i in range(6)])
            e = hs(v,c["sing"]*v).real/6
            if f==0 or sign==-1:
                assert e < cut
                if f == -1:
                    continue
                out.append([v.real for v in z])
                if f:
                    out.append([v.imag for v in z])
            else:
                assert e > cut
    assert len(out)==4
    return out


def fields(cut,n):
    zero = acb(0)
    out = [[[s[i]]+[zero]*8 for i in range(6)] for s in singlet_band(cut)]
    fam = real_family_kernel(cut,n)
    for j in range(fam.ncols()):
        out.append([[zero]+[fam[8*i+b,j] for b in range(8)] for i in range(6)])
    return out


def output_basis():
    l,r = (m["yuk"]["zero_basis"](sign) for sign in (1,-1))
    gen = [m["yuk"]["identity"](3)/arb(3).sqrt()]+c["gen"]
    out = []
    for i in range(6):
        L = acb_mat([[l[3*i+j,b] for b in range(3)] for j in range(3)])
        R = acb_mat([[r[3*i+j,b] for b in range(3)] for j in range(3)])
        out.extend([p[i]*c["phases"][i]*L.transpose()*g.transpose()*R for g in gen])
    return acb_mat([v.entries() for v in out]).transpose()


if __name__ == "__main__":
    lo,hi = arb("4.7983523138"),arb("4.7983523141")
    assert [count(rows(c["family"]),x) for x in (lo,hi)] == [30,31]
    assert [count(rows(c["sing"]),x) for x in (lo,hi)] == [4,4]
    assert len(fields(lo,30))==16
    vectors = fields(hi,31)
    assert len(vectors)==18
    X = acb_mat([[v for row in f for v in row] for f in vectors]).transpose()
    Y = output_basis()*X
    paired = acb_mat(rows(Y)+rows(Y.conjugate()))
    G,S = adj(paired)*paired,adj(X)*X
    assert count(rows(G-arb(".000119")**2*S),arb(0))==0
    print("first universal paired cutoff",lo.str(14),hi.str(14),
          "simultaneous gain > .000119; common complex dimensions 16 -> 18")
