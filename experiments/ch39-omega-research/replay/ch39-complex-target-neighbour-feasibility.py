"""中央复味空间的两个物理目标在相邻实际谱投影上的有界可达性。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, arb_mat


def setup(root):
    center = 2_306_888_100_085
    case = run_path(str(root/"ch39-complex-four-site-unit-quark-amplitudes.py"))
    native = run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    frame, weights = case["verify"](root)
    rows = {k:native["native"](root,center+k) for k in (-1,0,1)}
    targets = [sum((w*p for p,w in zip(rr,ws,strict=True)),acb_mat(4,4))
               for rr,ws in zip(rows[0],weights,strict=True)]
    return case,frame,rows,targets


def injection(frame, adj):
    basis = []
    for j in range(6):
        p,r,x,y,w,v = [int(j==k) for k in range(6)]
        z = acb(x,y)
        h = acb_mat([[p,z,0,0],[z.conjugate(),r,0,0],
                     [0,0,w,acb(0,v)],[0,0,acb(0,-v),w]])
        basis.append(adj(frame)*h*frame)
    gram = arb_mat([[(a*b).trace().real for a in basis] for b in basis])
    determinant = gram.det()
    assert determinant > 0
    print("固定复标架的六方向 Gram 行列式",determinant.str(12))


def separate(rows, target, offsets, coefficients, coord):
    def evaluate(p, vector):
        trace_w = (p[2,2]+p[3,3]).real
        return sum(x*z for x,z in zip(coord(p),vector,strict=True))-arb("1e-5")*trace_w

    values = [sum(evaluate(sector[j],form) for sector,form in zip(rows,coefficients,strict=True))
              for j in range(19)]
    goal = sum(evaluate(target,form) for form in coefficients)
    maximum = max(x.upper() for x in values)
    assert maximum < 0 < goal
    print("偏移、整数泛函、列值上界、目标值",offsets,coefficients,maximum.str(13),goal.str(13))


def verify(root):
    case,frame,rows,targets=setup(root)
    adj=case["adj"]
    injection(frame,adj)
    rotation=case["four_basis"]()
    right=rows[1][1]
    target=targets[1]
    _,weights=case["reconstruct"](right,adj(rotation)*target*rotation,rotation,[5,6,9,11])
    raw=sum((w*p for w,p in zip(weights,right,strict=True)),acb_mat(4,4))
    assert (adj(frame)*(raw-target)*frame).contains(acb_mat(3,3))
    assert all(0<w<arb(".396") for w in weights)
    print("右邻下型内点、按能序十九个权",*[v.str(11) for v in weights])
    for offsets,sector,vector in (
        ((-1,),0,[(-10,0,2,3)]),
        ((-1,),1,[(-91,7,20,-35)]),
        ((1,),0,[(-272,19,125,-26)]),
        ((0,1),1,[(-69,9,-30,8),(61,-5,23,-1)]),
    ):
        separate([rows[i][sector] for i in offsets],targets[sector],offsets,vector,case["coordinate_t"])


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
