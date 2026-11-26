"""Restore cyclic relabeling symmetry with three branches, using the original matrix field."""
from itertools import combinations
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, arb_mat

six = run_path(str(Path(__file__).with_name("ch39-shell-gw-six-well-selection.py")))
m = six["m"]
hs,adj,count,rows = (six[x] for x in ("hs","adj","count","rows"))
q = six["q"]
model = m["m"]["c"]["c"]
basis = [m["m"]["m"]["eye"](3)/arb(3).sqrt()]+model["gen"]
P = acb_mat([[int(i==(j+1)%3) for j in range(3)] for i in range(3)])


def certify(z):
    shapes = [sum((z[9*i+b,0]/q[i]*a for b,a in enumerate(basis)),acb_mat(3,3)) for i in range(6)]
    orbit = shapes+[P*s*adj(P) for s in shapes]+[adj(P)*s*P for s in shapes]
    differences = [v-orbit[-1] for v in orbit[:-1]]
    G = arb_mat([[hs(a,b).real for b in differences] for a in differences])
    # Change coordinates before the inertia test to keep interval errors small.
    numeric = np.array([[float(G[i,j]) for j in range(17)] for i in range(17)])
    _,vectors = np.linalg.eigh(numeric)
    V = arb_mat(vectors.T.tolist())
    metric = V*V.transpose()
    assert count(rows(acb_mat(metric)),arb(".99"))==0
    test = V*G*V.transpose()
    assert count(rows(acb_mat(test-arb(".0000129")*metric)),arb(0))==0
    omega = acb(0,2*arb.pi()/3).exp()
    target = [omega**(j//6)-omega**2 for j in range(17)]
    data = arb_mat([[z.real,z.imag] for z in target])
    transformed = V*data
    cov = transformed.transpose()*test.solve(transformed)
    gain = cov[0,0].sqrt()
    vertices = set(range(6))
    edges = [(i,j) for i,j in combinations(vertices,2) if i!=(j^1)]
    boundary = min(sum((i in part)!=(j in part) for i,j in edges)
                   for size in range(1,6) for part in combinations(vertices,size))
    assert boundary==4
    energy = hs(z,m["K"]*z).real+six["lam"]*m["energy_ref"](q)+sum(six["depths"])
    gap = 6*min(six["depths"])+3*boundary*26-energy
    assert gap > arb("6.38")
    print("18 affine wells, real Gram > .0000129, orbit coordinate norm",gain.str(8),
          "mixed cyclic branches exceed ground energy by",gap.lower().str(9))


if __name__=="__main__":
    z = m["exact_fit"](acb_mat([[acb(x,y)] for x,y in six["seed"]]),six["t"],six["lam"])
    certify(z)
