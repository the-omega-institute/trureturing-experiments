"""三枚相邻原壳零点、两种夸克荷的同一实谱多项式插值秩。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, arb_mat, ctx


def moment_map(matrix,positions):
    # 在四站循环不变压缩中读取六个实坐标：极点对角、射线对角、两类复边。
    seed=acb_mat([[int(i==j) for j in positions] for i in range(19)])
    current=seed
    center=acb_mat([[int(i==j)*17 for i in range(19)] for j in range(19)])
    step=(matrix-center)*arb(1)/16
    cols=[]
    for _ in range(36):
        pole,first=positions[:2]
        rayedge=current[first,2]
        spoke=current[pole,1]
        cols.append([current[pole,0].real,current[first,1].real,spoke.real,spoke.imag,
                     rayedge.real,rayedge.imag])
        current=step*current
    return arb_mat(cols).transpose().tolist()


def family_rank(root):
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    flavor=run_path(str(root/"ch39-three-gram-source-flavor.py"))
    frame,_=flavor["family_space"](native["gram_cholesky"](root))
    adj=flavor["adj"]
    generators=[]
    for j in range(6):
        p,r,s,t,u,v=[int(i==j) for i in range(6)]
        z,w=acb(s,t),acb(u,v)
        m=acb_mat([[p,z,z,z],[z.conjugate(),r,w,w.conjugate()],
                   [z.conjugate(),w.conjugate(),r,w],[z.conjugate(),w,w.conjugate(),r]])
        generators.append(adj(frame)*m*frame)
    det=arb_mat([[(a*b).trace().real for b in generators] for a in generators]).det()
    print("三源上六个厄米方向的 Gram 行列式",det.str(13),flush=True)
    assert arb("148.99478")<det<arb("148.99480")


def verify(root):
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec=416
    pos=[shell["soft"]["pole"],*shell["boundary_sites"]()]
    assert pos==[14,8,7,6]
    table=[]
    n=2_306_888_100_085
    for j in (n-1,n,n+1):
        height=acb.zeta_zero(j).imag
        for charge in (-4,2):
            lap=shell["laplacian"](charge,height)
            table.extend(moment_map(lap,pos))
    matrix=arb_mat(table)
    midpoint=arb_mat([[z.mid() for z in row] for row in table])
    mag=midpoint.det()
    print("中点行列式",mag.str(8),flush=True)
    inv=midpoint.inv()
    err=inv*matrix-arb_mat([[int(i==j) for j in range(36)] for i in range(36)])
    norm=max(sum(abs(z).upper() for z in row) for row in err.tolist())
    print("三十六列预条件缺陷",norm.str(9),flush=True)
    assert norm<arb(".032")
    family_rank(root)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
