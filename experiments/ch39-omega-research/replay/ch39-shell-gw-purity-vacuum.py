"""Global singlet bounds and the local loss of stability for a fixed-length matrix Higgs."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

m = run_path(str(Path(__file__).with_name("ch39-shell-gw-higgs-full-band.py")))
s, k = m["sing"], m["family"]
hs, rows, count = (m[x] for x in ("hs", "rows", "count"))


def ground():
    v = m["frame"](s, arb(".02"))[0]
    return (hs(v, s*v)/hs(v, v)).real, [v[i,0]/abs(v[i,0]) for i in range(v.nrows())]


def schur_matrix():
    lam, phases = ground()
    assert arb(".01375939") < lam < arb(".01375941")
    assert [count(rows(k), arb(x)) for x in (".01588062", ".01588064", ".01588312")] == [0, 1, 2]
    assert [count(rows(s), arb(x)) for x in (".01375939", ".01375941")] == [0, 1]
    re = acb_mat(k.nrows(), k.ncols())
    im = acb_mat(k.nrows(), k.ncols())
    for i in range(k.nrows()):
        for j in range(k.ncols()):
            z = phases[i//8].conjugate()*k[i,j]*phases[j//8] - (lam if i==j else 0)
            re[i,j], im[i,j] = z.real, z.imag
    assert count(rows(re), arb("0")) == 0
    imag_lift = -re.solve(im)
    return (re-im*imag_lift)*arb(3)/4, imag_lift, lam, phases


def array(mat):
    return np.array([[complex(x).real for x in row] for row in rows(mat)])


def finite_witness(effective, lift, phases, lam):
    _, vectors = np.linalg.eigh(array(effective))
    x = vectors[:, :2]
    modes = (x + 1j * array(lift) @ x).T.reshape(2,6,8)
    generators = m["g"]
    eye = m["yuk"]["identity"](3)
    # Pick the angle by the size and sign of the leading cubic purity gain.
    gen = np.array([[[complex(z) for z in row] for row in rows(g)] for g in generators])
    def field(theta):
        z = np.cos(theta)*modes[0]+np.sin(theta)*modes[1]
        return z * (.333333 / max(np.linalg.norm(z,axis=1)))
    def cubic(z):
        h = np.einsum("va,aij->vij", z, gen)
        ad = h.conj().transpose(0,2,1)
        return np.trace((h+ad)@ad@h, axis1=1, axis2=2).real.sum()
    trial = max((field(theta) for theta in np.linspace(0,2*np.pi,361)), key=cubic)
    fields = []
    scalars = []
    purity = arb(0)
    for i, row in enumerate(trial):
        coords = [acb(float(z.real),float(z.imag)) for z in row]
        magnitude = (1-sum(abs(z)**2 for z in coords)).sqrt()
        scalar = phases[i]*magnitude
        scalars.append([scalar])
        fields.extend([[phases[i]*c] for c in coords])
        h = magnitude*eye/arb(3).sqrt() + sum((c*g for c,g in zip(coords,generators)),acb_mat(3,3))
        traceless = m["adj"](h)*h - eye/3
        purity += hs(traceless,traceless).real
    a, b = acb_mat(scalars), acb_mat(fields)
    kinetic = hs(a,s*a).real+hs(b,k*b).real-6*lam
    return kinetic/purity


if __name__ == "__main__":
    effective, lift, lam, phases = schur_matrix()
    for threshold, n in ((".0015918842",0),(".0015918845",2)):
        assert count(rows(effective),arb(threshold)) == n
    global_lower = (arb(".01588062")-lam)*2/3
    global_upper = finite_witness(effective,lift,phases,lam)
    assert global_lower > arb(".00141414")
    assert global_upper < arb(".0014151")
    print("global crossing bounded by",global_lower.lower().str(9),
          global_upper.upper().str(9),"local loss of stability in (.0015918842,.0015918845)")
