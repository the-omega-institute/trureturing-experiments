"""第三步只转动七维平凡循环空间，在 N+4 精确实现三个上型质量中心。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    full = lib("next-common-ckm-gravity")
    cx,rx,cy,ry,_ = full["verify"](root)
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    emb = lib("common-embedding-central-right-ckm")
    _,labels,_ = emb["verify_nodes"](root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    source,target = [spec(shell,-4,node)[1] for node in (n+1,n+4)]
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"]
    x0 = seed(root,rank["eye"],"ch39-full-ckm-gravity-cayley.tsv")
    y0 = seed(root,rank["eye"],"ch39-next-common-ckm-gravity-cayley.tsv")
    w = labels[n+1,-4]
    ov = adj(target)*adj(y0*y0)*adj(x0*x0)*source
    density = acb_mat([[ov[j,i]*w[i] for i in range(19)] for j in range(19)])*adj(ov)
    geo = shell["geo"]["sites"]
    rotate = [geo.index((y,z,x)) for x,y,z in geo]
    assert (19+2*sum(i==j for i,j in enumerate(rotate)))//3==7
    trivial = [0,3,6,9,11,15,17]
    assert rotate[14]==14
    for j in trivial:
        assert abs(target[14,j])>0
        assert all((target[i,j]-target[rotate[i],j]).contains(0) for i in range(19))
    sites = (shell["soft"]["pole"],*shell["boundary_sites"]())
    assert sites==(14,8,7,6)
    four = acb_mat([[target[i,j] for j in range(19)] for i in sites])
    sub = acb_mat([[density[i,j] for j in trivial] for i in trivial])
    eig = sub.eig(algorithm="rump")
    assert all(z.imag.contains(0) for z in eig)
    eig = sorted((z.real for z in eig),key=lambda v:v.mid())
    assert eig[0]>0 and all(a<b for a,b in zip(eig,eig[1:]))
    # 平凡块按序排列后，仅混合最高两个源本征值。其余十二维上取恒等。
    average = (eig[5]+eig[6])/2
    assigned = (eig[2],average,average,eig[3],eig[4],eig[0],eig[1])
    weights = [density[i,i].real for i in range(19)]
    for i,weight in zip(trivial,assigned,strict=True):
        weights[i] = weight
    raw = acb_mat([[four[i,j]*weights[j] for j in range(19)] for i in range(4)])*adj(four)
    vals,_ = lib("full-ckm-generator-next-shell")["sorted_eigen"](raw)
    def seed_radius(center,rad,directions):
        assert len(set(directions))==len(directions)==len(center)
        sizes = [abs(z+arb(0,rad.upper())).upper() for z in center]
        return 2*sum(p*p*(1 if i==j else 2) for p,(i,j,_) in zip(sizes,directions,strict=True)).sqrt()
    dx = seed_radius(cx,rx,lib("one-generator-gravity-full-ckm")["DIRECTIONS"])
    dy = seed_radius(cy,ry,full["DIRECTIONS"])
    scale = max(z.upper() for z in w)-min(z.lower() for z in w)
    error = 2*(dx+dy)*scale
    print("有理种子的七维压缩谱",[z.str(12) for z in eig],flush=True)
    print("种子四站谱、覆盖两个参数盒的统一误差上界",[z.str(14) for z in vals],error.upper().str(12),flush=True)
    wiggle = arb(0,error.upper())
    levels = [v+wiggle for v in vals]
    mass = list(map(arb,rank["CENTERS"]))
    assert all(levels[i]<m<levels[i+1] for i,m in enumerate(mass))
    prob = [prod(v-m for m in mass)/prod(v-z for j,z in enumerate(levels) if i!=j) for i,v in enumerate(levels)]
    assert all(p>0 for p in prob) and sum(prob).contains(1)
    print("三中心相邻两侧的最小严格距离",min(min((m-levels[i]).lower(),(levels[i+1]-m).lower())
          for i,m in enumerate(mass)).str(14),flush=True)
    print("四站三平面法向概率包围",[p.str(10) for p in prob],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
