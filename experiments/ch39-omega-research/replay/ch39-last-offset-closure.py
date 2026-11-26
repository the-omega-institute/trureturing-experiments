"""Uniform closure using signed interpolation rules for one moving seventh offset."""
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
old = run_path(str(root/'ch39-eight-time-pole-control.py'))
B, node, rounded, mul = (old[k] for k in ('B','node','rounded','mul'))
tolerance = Q(1, 2**272)


def interpolate(ids, j, *, last_cell=None, error=tolerance):
    offsets = [1,2,3,4,5,6,j]
    squares = [k*k for k in offsets]
    nodes = [node(i) for i in ids]
    weights = []
    for i, v in enumerate(nodes):
        others = nodes[:i]+nodes[i+1:]
        coefficient = B(Q(1,prod(squares)))
        for k,r in zip(offsets,squares,strict=True):
            other_factors = prod(1-Q(s,r) for s in squares if r != s)
            residue = mul(*(B(Q(-1,r))-x for x in others))/B(other_factors)
            if k == j and last_cell is not None:
                t = last_cell
            else:
                v0 = Q(old['levels'][k],2**288)
                t = B(v0,v0+Q(1,2**288))
            coefficient += residue*B(t.lo-error,t.hi+error)
        weights.append(rounded(coefficient*mul(*(B(1)+B(r)*v for r in squares))
                               /mul(*(v-x for x in others))))
    return nodes,weights


def main():
    paired = [1,2,3,4,8,9,51,52]
    bounded = [0,1,2,3,5,6,15,16]
    assert all(j==i+1 for i,j in zip(paired[::2],paired[1::2]))
    assert all(j==i+1 for i,j in zip(bounded[4::2],bounded[5::2]))
    for j in range(7,17):
        a = interpolate(paired,j)
        d = interpolate(bounded,j)
        assert min(a[1][0].lo,d[1][2].lo,a[1][2].lo)>Q('.05')
        scale = rounded(a[1][0]*a[0][0]*old['b41'])
        ratio = old['boundary'](scale,(a,d,a),True)
        assert ratio.lo > Q('1.002')
        print(j,'lower ratio',float(ratio.lo))


if __name__ == '__main__':
    main()
