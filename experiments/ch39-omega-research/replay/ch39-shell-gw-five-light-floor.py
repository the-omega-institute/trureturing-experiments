"""Quantify how the two lighter singular values open the five-mode mixing gap."""
from itertools import combinations, combinations_with_replacement, permutations
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import acb_mat, arb

rank = run_path(str(Path(__file__).with_name("ch39-shell-gw-rank-one-mixing.py")))
model = rank["model"]
adj = rank["adj"]
pairs = list(combinations(range(3), 2))
quadratics = list(combinations_with_replacement(range(3), 2))
cubics = list(combinations_with_replacement(range(3), 3))


def minor_coefficients(comp):
    matrix = acb_mat(10, 4)
    for col, rows in enumerate(combinations(range(4), 3)):
        for perm in permutations(range(3)):
            sign = (-1)**sum(perm[i]>perm[j] for i,j in pairs)
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        idx = cubics.index(tuple(sorted((a,b,c))))
                        matrix[idx,col] += sign*comp[rows[0]][a,perm[0]]*comp[rows[1]][b,perm[1]]*comp[rows[2]][c,perm[2]]
    return matrix


def inverse_bounds(roots, comp):
    values = acb_mat([[x[i,0]*x[j,0] for x in roots] for i,j in quadratics])
    interpolation = values.inv()
    minor = minor_coefficients(comp)
    dual = (adj(minor)*minor).solve(adj(minor))
    bounds = []
    for i, x in enumerate(roots):
        product = acb_mat(10,3)
        for k,(a,b) in enumerate(pairs):
            for s,monomial in enumerate(quadratics):
                coeff = interpolation[i,s]
                ia = cubics.index(tuple(sorted((*monomial,a))))
                ib = cubics.index(tuple(sorted((*monomial,b))))
                product[ia,k] += x[b,0]*coeff
                product[ib,k] -= x[a,0]*coeff
        C = dual*product
        # The transpose acts on the vector of the four maximal minors.
        proposal = arb(float(np.linalg.svd(model["array"](C),compute_uv=False)[0]) + .00001)
        gram = adj(C)*C
        assert rank["seven"]["count"](rank["seven"]["rows"](gram),proposal**2)==3
        bounds.append(proposal)
    return bounds


def check(index):
    modes, _ = rank["seven"]["spectral_band"](index)
    left, right = rank["mix"]["orbit_frame"](modes)
    modes = [adj(left)*b*right for b in modes]
    comp = model["complement"](modes[:5])
    roots = [rank["certified_root"](comp,p,seed)[0] for p in range(3)
             for seed in rank["numeric_seeds"](comp,p)]
    norms = inverse_bounds(roots,comp)
    print(index,"whole-space bound",sum(norms).upper().str(7),
          "point bounds",[t.upper().str(6) for t in norms])
    distances = [[0 if i==j else (1-abs(model["hs"](v,w))**2).sqrt()
                  for j,w in enumerate(roots)] for i,v in enumerate(roots)]
    def bound(s):
        c = sum(norms)
        for _ in range(10):
            radius=c*s
            updated=[]
            for i,k in enumerate(norms):
                if any(not (d>radius) for j,d in enumerate(distances[i]) if j!=i):
                    return c
                loss=sum(b*s/(distances[i][j]-radius) for j,b in enumerate(norms) if j!=i)
                if not (loss<1):
                    return c
                updated.append(k/(1-loss))
            c=min(c,max(updated))
        return c
    min_dist=min(d for i,row in enumerate(distances) for j,d in enumerate(row) if i!=j)
    cap=arb(".0055" if index==43 else ".0049")
    factor=bound(cap)
    assert factor*cap < arb(".0428") and min_dist-factor*cap>arb(".043")
    print("certified light mass floor",cap.str(7),"angle sensitivity",factor.upper().str(7))
    return cap


if __name__ == "__main__":
    for n in (43,195):
        check(n)
