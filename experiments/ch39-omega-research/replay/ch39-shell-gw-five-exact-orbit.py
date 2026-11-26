"""Interval certificate that the five-mode space meets the exact reference CKM orbit."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, arb_mat

folder = Path(__file__).resolve().parent
four = run_path(str(folder/"ch39-shell-gw-four-mode-obstruction.py"))
five = run_path(str(folder/"ch39-shell-gw-higgs-five-fit.py"))
hs, adj = five["hs"], five["adj"]
as_array, exact_frame = four["float_matrix"], four["exact_frame"]
orthonormalize = four["orthonormalize"]
identity = five["yuk"]["identity"]


def reference():
    v = identity(3)
    for i, j, angle, ph in ((1,2,".04275","0"), (0,2,".00377","-1.139"), (0,1,".22704","0")):
        t = arb(angle)
        b = identity(3)
        b[i,i] = b[j,j] = t.cos()
        b[i,j] = acb(0,arb(ph)).exp()*t.sin()
        b[j,i] = -b[i,j].conjugate()
        v = v*b
    def diagonal(values):
        return acb_mat([[arb(values[i])/arb(values[2]) if i==j else 0
                         for j in range(3)] for i in range(3)])
    return v, diagonal((".00000615",".00311",".8616")), diagonal((".0000135",".000268",".01401"))


def trial_coordinates(ckm):
    def svd(target):
        y = np.array([complex(*(float(c) for c in x.split())) for x in target]).reshape(3,3)
        l, _, rh = np.linalg.svd(y)
        return l[:, ::-1], rh[::-1, :]
    u, a = svd(five["up"])
    v, b = svd(five["down"])
    m = (u.conj().T @ v) / as_array(ckm)
    p = m[:,2]/abs(m[:,2])
    q = m[0,:]/p[0]
    q = q/abs(q)
    return tuple(exact_frame(x) for x in (u@np.diag(p), np.diag(p.conj())@a, np.diag(q)@b))


def complement():
    _, ys = five["matrices"]()
    modes = orthonormalize(ys)
    basis = np.stack([as_array(b).flatten() for b in modes]).T
    trials = np.linalg.svd(basis, full_matrices=True)[0][:,5:]
    vectors = []
    for v in trials.T:
        t = acb_mat([[acb(float(z.real), float(z.imag)) for z in row] for row in v.reshape(3,3)])
        for b in modes:
            t = t-hs(b,t)*b
        vectors.append(t)
    return orthonormalize(vectors)


def generators():
    out = []
    for i in range(3):
        d = acb_mat(3,3)
        d[i,i] = acb(0,1)
        out.append(d)
        for j in range(i):
            for z in (acb(1),acb(0,1)):
                m = acb_mat(3,3)
                m[i,j] = z/arb(2).sqrt()
                m[j,i] = -m[i,j].conjugate()
                out.append(m)
    return out


def real_frame_pair(comp, x, y):
    zs = [hs(c,m) for m in (x,y) for c in comp]
    return [z.real for z in zs]+[z.imag for z in zs]


def frob(a):
    return sum(a[i,j]**2 for i in range(a.nrows()) for j in range(a.ncols())).sqrt()


def exact_certificate():
    ck, du, dd = reference()
    common, ru, rd = trial_coordinates(ck)
    x = common*du*ru
    y = common*ck*dd*rd
    comp = complement()
    proj = lambda a,b: real_frame_pair(comp,a,b)
    residual = arb_mat([[v] for v in proj(x,y)])
    zero = acb_mat(3,3)
    gen = generators()
    cols = [proj(g*x,g*y) for g in gen]+[proj(x*g,zero) for g in gen]+[proj(zero,y*g) for g in gen]
    J = arb_mat([list(c) for c in zip(*cols)])
    # The exact right inverse C = Jᵀ (JJᵀ)⁻¹ parametrizes a transverse slice.
    C = J.transpose()*(J*J.transpose()).inv()
    n, eta = frob(C), frob(residual)
    assert n < 10000 and eta < arb("1e-10")
    radius = 2*eta
    # In right-inverse coordinates the derivative at the origin is the identity.
    # The unitaries are exponentials of skew-Hermitian matrices for every real z.
    # Duhamel's formula bounds both exponential derivatives without growth factors.
    hessian = 6*n**2
    assert hessian*radius < arb(".5")
    # Perturb the two fixed normalized input matrices in the product HS metric.
    eps = arb("1.5e-6")
    radius_perturbed = 2*(eta+eps)
    slope_error = 2*n*eps
    second = hessian*(1+eps)
    assert slope_error+second*radius_perturbed < 1
    assert (eta+eps)/radius_perturbed+slope_error+second*radius_perturbed/2 < 1
    print("exact central inputs; right-inverse norm below", n.upper().str(9),
          "initial residual below",eta.upper().str(9),
          "contraction bound", (hessian*radius).upper().str(7))
    return n, eta


if __name__ == "__main__":
    exact_certificate()
