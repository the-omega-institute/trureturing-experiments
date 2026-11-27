"""跨七维平凡空间的一次微小旋转在 N+4 同时实现三个下型中心。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(m):
    return m.conjugate().transpose()


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    full = lib("next-common-ckm-gravity")
    cx,rx,cy,ry,_ = full["verify"](root)
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    _,labels,_ = lib("common-embedding-central-right-ckm")["verify_nodes"](
        root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    eye = rank["eye"](19)
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"]
    local = lib("next-common-ckm")["local"]
    def at_center(name,center,radius,directions):
        m = local(seed(root,rank["eye"],name),center,directions,eye,False)[0]
        assert len(set(directions))==len(directions)==len(center)
        widths = [abs((z+arb(0,radius.upper()))-z).upper() for z in center]
        delta = 2*sum(p*p*(1 if i==j else 2) for p,(i,j,_) in zip(widths,directions,strict=True)).sqrt()
        return m,delta
    x,dx = at_center("ch39-full-ckm-gravity-cayley.tsv",cx,rx,lib("one-generator-gravity-full-ckm")["DIRECTIONS"])
    y,dy = at_center("ch39-next-common-ckm-gravity-cayley.tsv",cy,ry,full["DIRECTIONS"])
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    source,target = [spec(shell,2,v)[1] for v in (n+1,n+4)]
    w = labels[n+1,2]
    ov = adj(target)*y*x*source
    rho = acb_mat([[ov[i,j]*w[j] for j in range(19)] for i in range(19)])*adj(ov)
    error = (max(z.upper() for z in w)-min(z.lower() for z in w))*(dx+dy)
    wiggle = arb(0,error.upper())
    sites = shell["geo"]["sites"]
    cycle = [sites.index((b,c,a)) for a,b,c in sites]
    groups = lib("charge-six-projector-algebra-transport")["groups"](target,cycle)
    seven,other = groups[0],sorted(groups[1]+groups[2])
    def sub(rows,cols):
        return acb_mat([[rho[i,j] for j in cols] for i in rows])
    eig = lib("full-ckm-generator-next-shell")["sorted_eigen"]
    a,avec = eig(sub(seven,seven))
    b,bvec = eig(sub(other,other))
    cross = abs((adj(avec[4])*sub(seven,other)*bvec[5])[0,0])
    assert arb(".01770")<cross<arb(".01772")
    root3 = arb(3).sqrt()
    omega = acb(-arb(1)/2,root3/2)
    support = []
    for phase in (omega,omega.conjugate()):
        strengths = [abs(target[8,j]+phase*target[7,j]+phase**2*target[6,j])**2/3 for j in range(19)]
        used = sorted((j for j in other if strengths[j]>0),key=lambda j:strengths[j].mid(),reverse=True)
        assert len(used)==6 and all(strengths[i]>strengths[j] for i,j in zip(used,used[1:]))
        support.append(used)
    assert len(set(seven+support[0]+support[1]))==19
    print("平凡、第一字符和第二字符的目的层编号",[[j+1 for j in g] for g in (seven,*support)],flush=True)
    keep = (shell["soft"]["pole"],*shell["boundary_sites"]())
    assert keep==(14,8,7,6)
    four = acb_mat([[target[i,j] for j in range(19)] for i in keep])
    masses = tuple(map(arb,(".0000135",".000268",".01401")))
    for angle in (".0066",".0063765"):
        sn,cs = arb(angle).sin_cos()
        moved_a = a[4]*cs**2+b[5]*sn**2-2*cross*sn*cs
        moved_b = a[4]+b[5]-moved_a
        va,vb = a.copy(),b.copy()
        va[4],vb[5] = moved_a,moved_b
        masses7 = [va[k] for k in (2,5,3,4,6,1,0)]
        values = [arb(0)]*19
        for group,weight in zip((seven,*support),(masses7,vb[:6],vb[6:]),strict=True):
            for j,val in zip(group,weight,strict=True):
                values[j]=val
        raw = acb_mat([[four[i,j]*values[j] for j in range(19)] for i in range(4)])*adj(four)
        levels = [z+wiggle for z in eig(raw)[0]]
        assert all(levels[i]<m<levels[i+1] for i,m in enumerate(masses))
        p = [prod(z-m for m in masses)/prod(z-y for j,y in enumerate(levels) if j!=i)
             for i,z in enumerate(levels)]
        assert all(z>0 for z in p) and sum(p).contains(1)
        print("角度、泄漏概率、跨块矩阵元及整盒谱误差上界",angle,(sn*sn).str(12),
              cross.str(12),error.upper().str(12),flush=True)
        print("四站有序谱和法向概率",[z.str(13) for z in levels],[z.str(10) for z in p],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
