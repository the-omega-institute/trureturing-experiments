"""中心十一阶共同味多项式冻结后，在两枚相邻零点检验平方质量。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, arb_mat, ctx


def family_coordinates(frame,original):
    cols=[]
    for index in range(6):
        p,r,s,t,u,v=[int(index==j) for j in range(6)]
        z,w=acb(s,t),acb(u,v)
        block=acb_mat([[p,z,z,z],[z.conjugate(),r,w,w.conjugate()],
                       [z.conjugate(),w.conjugate(),r,w],[z.conjugate(),w,w.conjugate(),r]])
        cols.append(original["coordinates"](original["adj"](frame)*block*frame))
    return arb_mat(cols).transpose()


def moments(ray, positions, number):
    height=acb.zeta_zero(number).imag
    identity=acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    seeds=acb_mat([[int(i==j) for j in positions] for i in range(19)])
    rows=[]
    for charge in (-4,2):
        lap=ray["laplacian"](charge,height)
        scaled=(lap-18*identity)*(arb(1)/8)
        current=seeds
        columns=[]
        for _ in range(12):
            p,r=positions[:2]
            z,w=current[p,1],current[r,2]
            columns.append([current[p,0].real,current[r,1].real,z.real,z.imag,w.real,w.imag])
            current=scaled*current
        rows.extend(arb_mat(columns).transpose().tolist())
    return arb_mat(rows)


def verify(root):
    original=run_path(str(root/"ch39-signed-three-source-quark-pair.py"))
    ray=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    ckm=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    flavor=run_path(str(root/"ch39-three-gram-source-flavor.py"))
    ctx.prec=320
    frame=original["adapted_frame"](native["gram_cholesky"](root))
    family=family_coordinates(frame,original)
    targets=original["target_blocks"](ckm,flavor)
    raw=[]
    for block in targets:
        target=arb_mat([[v] for v in original["coordinates"](block)])
        raw.extend(family.solve(target).tolist())
    rhs=arb_mat(raw)
    sites=[ray["soft"]["pole"],*ray["boundary_sites"]()]
    assert sites==[14,8,7,6]
    n=2_306_888_100_085
    central=moments(ray,sites,n)
    assert not central.det().contains(0)
    coefficients=central.solve(rhs,algorithm="precond")
    minima={(n-1,-4):("-0.505","-0.504"),(n-1,2):("-601.825","-601.824"),
            (n+1,-4):("-3550.395","-3550.394"),(n+1,2):("-1.992","-1.991")}
    for surrounding in (n-1,n+1):
        other=moments(ray,sites,surrounding)
        direct=other*coefficients
        dual=central.transpose().solve(other.transpose(),algorithm="precond").transpose()*rhs
        joint=[a.intersection(b) for a,b in zip(direct.entries(),dual.entries(),strict=True)]
        for j,charge in enumerate((-4,2)):
            coordinate=family*arb_mat([[v] for v in joint[6*j:6*j+6]])
            a,b,c,s,t,h=coordinate.entries()
            z,w=acb(s,t),acb(0,h)
            matrix=acb_mat([[a,0,z],[0,b,w],[z.conjugate(),w.conjugate(),c]])
            intervals=matrix.eig(algorithm="rump")
            assert all(v.imag.contains(0) for v in intervals)
            levels=sorted((v.real for v in intervals),key=lambda v:float(v.mid()))
            assert levels[0]<levels[1]<levels[2]
            lo,hi=map(arb,minima[surrounding,charge])
            assert lo<levels[0]<hi and levels[2]>0
            expected=1 if (surrounding,charge)==(n-1,-4) else 2
            assert sum(v<0 for v in levels)==expected and all(not v.contains(0) for v in levels)
            print("留出编号、荷、三个平方质量",surrounding,charge,*[v.str(12) for v in levels],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
