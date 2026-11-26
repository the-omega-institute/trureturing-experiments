"""Scalar and gauge metric tests on periodic curved simplex geometries."""
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, fmpq, fmpq_mat
import sympy as sp

base = run_path(str(Path(__file__).with_name("ch39-scalar-block-metric.py")))
planes = tuple(combinations(range(4),2))


def wedge(M):
    return type(M)([[M[i,k]*M[j,l]-M[i,l]*M[j,k] for k,l in planes] for i,j in planes])


def circulation(order,parity,parities):
    masks = [sum(1<<v for v in order[:i]) for i in range(5)]
    P = fmpq_mat([[(m>>i)&1 for i in range(4)] for m in masks[1:]])
    edges = 15*len(parities)
    rows = []
    for a,b in planes:
        row = [0]*edges
        for i,j,s in [(0,a+1,1),(a+1,b+1,1),(0,b+1,-1)]:
            mask = masks[j]-masks[i]
            start_parity = parity*(-1)**i
            row[15*parities.index(start_parity)+mask-1] += 2*s
        rows.append(row)
    return wedge(P).inv()*fmpq_mat(rows)


def submatrix(M,rows,cols):
    return type(M)([[M[i,j] for j in cols] for i in rows])


def stiffness(parities):
    mats,rows = [],[]
    scalar = arb_mat(4,4)
    for p in parities:
        for order in base["orders"]:
            det,_,inv = base["simplex"](order,p)
            vol = arb(det).sqrt()/(24*len(parities))
            C = circulation(order,p,parities)
            rows.extend(C.tolist())
            scalar += arb_mat(inv)*vol
            mats.append((vol,wedge(inv),C))
    reduced,rank = fmpq_mat(rows).rref()
    pivots = [next(j for j in range(reduced.ncols()) if reduced[i,j]) for i in range(rank)]
    assert rank == (25 if len(parities)==2 else 11)
    H = arb_mat(6+rank,6+rank)
    for vol,W,C in mats:
        E = fmpq_mat([[int(i==j) for j in range(6)]+[C[i,j] for j in pivots] for i in range(6)])
        H += arb_mat(E.transpose()*W*E)*vol
    field,adjust = range(6),range(6,6+rank)
    cross = submatrix(H,field,adjust)
    cell = submatrix(H,field,field)-cross*submatrix(H,adjust,adjust).inv()*cross.transpose()
    return wedge(scalar)/scalar.det().sqrt(),cell


def test_ratio(K,A,B):
    def energy(v):
        vector = arb_mat([v])
        return (vector*K*vector.transpose())[0,0]
    return energy(A)/energy(B)


def verify():
    # Use the same two constant two-forms for all three probes.
    A = [1,0,0,0,0,0]
    B = [0,1,-1,-1,1,0]  # (e1-e2) wedge (e3-e4)
    equal = wedge(arb_mat([[1+int(i==j) for j in range(4)] for i in range(4)]))
    bounds = [[(74767,74768),(74853,74854)],[(68969,68970),(68947,68948)]]
    for parities,limits in zip(((-1,1),(0,)),bounds):
        scalar,whitney = stiffness(parities)
        ratios = [test_ratio(K,A,B) for K in (scalar,equal,whitney)]
        assert ratios[1] == fmpq(3,4)
        for value,(lo,hi) in zip((ratios[0],ratios[2]),limits):
            assert fmpq(lo,100000) < value and value < fmpq(hi,100000)
        print(parities,"scalar, equal faces, relaxed Whitney energy ratios:",*ratios)


def rational_directional_ratios():
    # Only the e1 edges change. Every simplex perturbation has rank at most two.
    xs = (-fmpq(1,3),fmpq(0),fmpq(1,3))
    vecs = ([1,0,0,0,0,0],[0,1,-1,-1,1,0])
    x = sp.symbols("u")
    sums = sp.zeros(4)
    gauge_matrix = sp.zeros(6)
    for order in base["orders"]:
        pos = order.index(0)
        P = fmpq_mat([[int(j in order[:i]) for j in range(4)] for i in range(1,5)])
        g0 = fmpq_mat([[int(i==j)-fmpq(1,5) for j in range(4)] for i in range(4)])
        grad = fmpq_mat([[-1]*4]+[[int(i==j) for j in range(4)] for i in range(4)])*P.inv().transpose()
        U = fmpq_mat([[ -(grad[pos,i]*grad[pos+1,j]+grad[pos+1,i]*grad[pos,j])/2
                        for j in range(4)] for i in range(4)])
        vals = []
        for v in xs:
            g = g0+U*v
            d = 1+v-3*v**2/4
            assert g.det() == d/5
            vals.append(sp.Matrix([[sp.Rational(str(z)) for z in row] for row in (g.inv()*d).tolist()]))
        neg,center,posval = vals
        inv = center+3*x*(posval-neg)/2+9*x*x*(posval+neg-2*center)/2
        sums += inv
        K = sp.Matrix([[inv.extract(a,b).det() for b in planes] for a in planes])
        gauge_matrix += K
    K = sp.Matrix([[sums.extract(a,b).det() for b in planes] for a in planes])
    gauge = [(sp.Matrix(v).T*gauge_matrix*sp.Matrix(v))[0] for v in vecs]
    scalar = [(sp.Matrix(v).T*K*sp.Matrix(v))[0] for v in vecs]
    r_s,r_w = sp.factor(scalar[0]/scalar[1]),sp.factor(gauge[0]/gauge[1])
    assert sp.cancel(r_s-sp.Rational(3,4)+5*x*(x-4)/(24*(x**2-x-4))) == 0
    factor = x*x*(x*x-18*x+8)
    assert sp.cancel(r_s-r_w - 5*factor/(24*(x*x-x-4)*(5*x*x-16*x-48))) == 0
    for i,j in combinations_with_replacement(range(6),2):
        cross = sp.Poly(K[i,j]*gauge_matrix[0,0]-K[0,0]*gauge_matrix[i,j],x)
        assert cross.rem(sp.Poly(factor,x)).is_zero
    special = 9-arb(73).sqrt()
    det = (1+special-3*special**2/4)/5
    assert det > (fmpq(4,5)+special)**4/16 and det > fmpq(6,5)**4/16
    print("Coordinate-edge pencil:", r_s, r_w, "difference", sp.factor(r_s-r_w))


if __name__ == "__main__":
    verify()
    rational_directional_ratios()
