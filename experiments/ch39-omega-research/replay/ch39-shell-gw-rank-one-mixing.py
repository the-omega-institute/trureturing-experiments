"""Recover all rank-one points of the five-mode space, including their principal angles."""
from itertools import combinations
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

mix = run_path(str(Path(__file__).with_name("ch39-shell-gw-three-six-mixing.py")))
model, seven, adj = mix["model"], mix["seven"], mix["adj"]
pairs = ((0, 1, 2), (0, 1, 3))


def sliced_minors(matrices):
    # Restrict x^T C to the projective line x=(z,w,z+iw). The two
    # selected maximal minors have no common projective zero on this line.
    cubics = []
    for triple in pairs:
        def value(z, w):
            x = acb_mat([[z, w, z+acb(0, 1)*w]])
            return acb_mat([(x*matrices[j]).entries() for j in triple]).det()
        a, d = value(1, 0), value(0, 1)
        plus, minus = value(1, 1), value(1, -1)
        cubics.append([a, (plus-minus)/2-d, (plus+minus)/2-a, d])
    sylvester = acb_mat([[0]*i + c + [0]*(2-i) for c in cubics for i in range(3)])
    assert abs(sylvester.det()) > 0


def numeric_seeds(matrices, pivot):
    # The two roots sought in this chart have the same dominant left/right axis.
    c = np.array(list(map(model["array"], matrices)))
    free = [j for j in range(3) if j != pivot]
    rng = np.random.default_rng(192)
    u = np.ones((32, 3), complex)
    v = u.copy()
    for vectors in (u, v):
        vectors[:, free] = rng.normal(size=(32, 2))+1j*rng.normal(size=(32, 2))
    for _ in range(40):
        residual = np.einsum("ti,gij,tj->tg", u, c, v)
        left = np.einsum("gij,tj->tgi", c[:, free, :], v)
        right = np.einsum("ti,gij->tgj", u, c[:, :, free])
        step = np.linalg.solve(np.concatenate([left,right],axis=2),residual[:,:,None])[:,:,0]
        u[:, free] -= step[:, :2]
        v[:, free] -= step[:, 2:]
    seeds = []
    for x,y in zip(u,v):
        if max(np.max(abs(x[free])), np.max(abs(y[free]))) >= .8:
            continue
        if all(np.linalg.norm(x-a) > .01 for a,b in seeds):
            seeds.append((x,y))
    assert len(seeds) == 2
    return seeds


def certified_root(matrices, pivot, seed):
    free = [j for j in range(3) if j != pivot]
    def axes(radius):
        error = arb(0,radius)
        return [acb_mat([[1 if i==pivot else acb(float(z.real),float(z.imag))+acb(error,error)]
                         for i,z in enumerate(v)]) for v in seed]
    def jacobian(u,v):
        return acb_mat([[(c*v)[j,0] for j in free]+[(u.transpose()*c)[0,j] for j in free]
                        for c in matrices])
    u,v = axes(0)
    res = acb_mat([[(u.transpose()*c*v)[0,0]] for c in matrices])
    inv = jacobian(u,v).inv()
    radius = arb("1e-10")
    assert max(abs(z) for z in (inv*res).entries()) < radius/2
    u,v = axes(radius)
    contraction = model["eye"](4)-inv*jacobian(u,v)
    assert max(sum(abs(contraction[i,j]) for j in range(4)) for i in range(4)) < arb(".1")
    return [x/model["hs"](x,x).real.sqrt() for x in (u,v)]


def check(index):
    modes, _ = seven["spectral_band"](index)
    left, right = mix["orbit_frame"](modes)
    modes = [adj(left)*b*right for b in modes]
    comp = model["complement"](modes[:5])
    sliced_minors(comp)
    sliced_minors([c.transpose() for c in comp])
    roots = [certified_root(comp,pivot,seed) for pivot in range(3)
             for seed in numeric_seeds(comp,pivot)]
    distances = []
    for (u,v),(x,y) in combinations(roots,2):
        distances.append((1-abs(model["hs"](u,x))**2).sqrt())
    lo,hi = (".3239",".3240") if index==43 else (".2930",".2931")
    assert all(d>arb(lo) for d in distances) and any(d<arb(hi) for d in distances)
    fifth = [abs((u.transpose()*modes[4]*v)[0,0]) for u,v in roots]
    assert all(a>arb("1e-6") for a in fifth)
    print("zero",index,"smallest nonzero five-mode projector distance",min(distances).str(9),
          "least necessary fifth coefficient",min(fifth).str(7))


if __name__=="__main__":
    for index in (43, 195):
        check(index)
