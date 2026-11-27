"""右侧原训练谱，无额外荷六约束时以同一酉发生元及固定三平面实现六个质量中心。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


DIRECTIONS = ((14, 4, 0), (17, 8, 1), (14, 3, 0), (13, 6, 0), (17, 7, 0), (12, 3, 0))
SEED = (
    ((425662017667, -378540469763), (421446923214, -384984486350), (422078344429, -379003190106)),
    ((710178156498, 34272146993), (-289939806912, 38215422258), (-288451711114, 35517117939)),
    ((-290856295966, 30465147743), (708972781523, 34414685553), (-289500085759, 31730913579)),
    ((-288779882732, 32741940088), (-288916234595, 36668241276), (712578587488, 33988290114)),
)


def adj(a):
    return a.conjugate().transpose()


def verify(root):
    rank = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    prior = run_path(str(root/"ch39-one-generator-gravity-up-strange-bottom.py"))
    train = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    spectrum = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    shells = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    n = rank["N"]+1
    _, labels, _ = train(root, ((n,".0294"),), ".1", (), ((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    eye = rank["eye"](19)
    center = prior["cayley"](root,rank["eye"],"ch39-six-centers-common-generator-cayley.tsv")
    seed = acb_mat([[acb(a,b) for a,b in row] for row in SEED])
    columns = []
    for j in range(3):
        v = acb_mat([[seed[i,j]] for i in range(4)])
        for old in columns:
            v -= old*(adj(old)*v)[0,0]
        norm = (adj(v)*v)[0,0].real.sqrt()
        assert norm > 0
        columns.append(v*(1/norm))
    frame = acb_mat([[v[i,0] for v in columns] for i in range(4)])
    assert (adj(frame)*frame).contains(rank["eye"](3))
    generators = []
    for i,j,imag in DIRECTIONS:
        v = acb_mat(19,19)
        v[i,j] = acb(0,1) if imag else acb(1)
        v[j,i] = -v[i,j].conjugate()
        generators.append(v)
    sectors = []
    for q in (-4,2):
        source,target = [spectrum(shells,q,number)[1] for number in (n,n+1)]
        short = acb_mat([[target[i,j] for j in range(19)] for i in (14,8,7,6)])
        sectors.append((q,source,target,adj(frame)*short,labels[n,q]))
    masses = [list(map(arb,row)) for row in (rank["CENTERS"],(".0000135",".000268",".01401"))]
    def system(params):
        k = sum((v*t for v,t in zip(generators,params,strict=True)),acb_mat(19,19))
        inv = (eye-k).inv()
        x = (2*inv-eye)*center
        partial = [2*inv*v*inv*center for v in generators]
        val,jac = [],[]
        for (q,source,target,c,weights),targets in zip(sectors,masses,strict=True):
            if q==-4:
                z = adj(x)
                u = z*z
                deriv = [adj(t)*z+z*adj(t) for t in partial]
            else:
                u,deriv=x,partial
            a = adj(target)*u*source
            def block(amps):
                return acb_mat([[c[i,j]*amps[j] for j in range(19)] for i in range(3)])*adj(c)
            matrix = block([sum(abs(a[j,k])**2*weights[k] for k in range(19)) for j in range(19)])
            directions = []
            for t in deriv:
                du = adj(target)*t*source
                directions.append(block([2*sum((a[j,k].conjugate()*du[j,k]).real*weights[k]
                                                for k in range(19)) for j in range(19)]))
            for m in targets:
                h = matrix-m*rank["eye"](3)
                # 三阶伴随矩阵以显式余子式求出，即使 m 正好是本征值也可计算。
                cofactor = acb_mat(3,3)
                for i in range(3):
                    for j in range(3):
                        sub = acb_mat([[h[k,l] for l in range(3) if l!=j] for k in range(3) if k!=i])
                        cofactor[i,j] = (-1)**(i+j)*sub.det()
                val.append(h.det().real)
                jac.append([sum((cofactor[i,j]*d[i,j]).real for i in range(3) for j in range(3)) for d in directions])
        return arb_mat([[v] for v in val]),arb_mat(jac)
    point=[arb(0)]*6
    for _ in range(2):
        val,jac=system(point)
        point=[(x-y).mid() for x,y in zip(point,jac.solve(val).entries(),strict=True)]
    radius=arb("1e-12")
    val,jac=system(point)
    pre=jac.mid().inv()
    _,wide=system([x+arb(0,radius.upper()) for x in point])
    defect=arb_mat([[int(i==j) for j in range(6)] for i in range(6)])-pre*wide
    contraction=max((abs(z)/radius+sum(abs(v) for v in row)).upper()
                    for z,row in zip((pre*val).entries(),defect.tolist(),strict=True))
    assert contraction<1,contraction
    print("六个质量方程的收缩比与参数盒半宽",contraction.str(12),radius.str(6))
    print("同一三平面上的六个精确质量",*[x.str(9) for row in masses for x in row])


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
