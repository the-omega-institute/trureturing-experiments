"""由两壳联合参数盒导出 N+4 下型七维与十二维压缩格点和四站射线强度下界。"""
import argparse
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


SCALE = 10**12
FILENAME = "ch39-fourth-shell-strange-cells.tsv"


def build(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    x,y,*_ = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    _,labels,_ = lib("common-embedding-central-right-ckm")["verify_nodes"](
        root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    source,raw_target = [spec(shell,2,k)[1] for k in (n+1,n+4)]
    adj = rank["adj"]
    # 固定每条谱线的整体相位，采用第一项模严格超过 1/10 的分量。
    phases = []
    for j in range(19):
        i = next(i for i in range(19) if abs(raw_target[i,j])>arb(".1"))
        assert all(abs(raw_target[k,j])<arb(".1") for k in range(i))
        pivot = raw_target[i,j]
        phases.append(pivot.conjugate()/abs(pivot))
    target = acb_mat([[raw_target[i,j]*phases[j] for j in range(19)] for i in range(19)])
    ov = adj(target)*y*x*source
    w = labels[n+1,2]
    assert min(w)>0
    H = acb_mat([[ov[j,i]*w[i] for i in range(19)] for j in range(19)])*adj(ov)
    sites = shell["geo"]["sites"]
    cycle = [sites.index((b,c,a)) for a,b,c in sites]
    groups = lib("charge-six-projector-algebra-transport")["groups"](target,cycle)
    blocks = (groups[0],sorted(groups[1]+groups[2]))
    rows = [("scale",SCALE),("entry_error",1),("shells",n+1,n+4),("sites",14,8,7,6),
            ("trial",1,0,-2,5),("trivial",*[i+1 for i in blocks[0]]),
            ("complement",*[i+1 for i in blocks[1]])]
    def grid(value):
        scaled = SCALE*value
        z = int((scaled.mid()+arb(1)/2).floor().unique_fmpz())
        assert abs(scaled-z)<1
        return z
    for indices in blocks:
        for i,src in enumerate(indices):
            for j in range(i,len(indices)):
                z = H[src,indices[j]]
                assert i!=j or z.imag.contains(0)
                rows.append(("cell",len(indices),i,j,grid(z.real),0 if i==j else grid(z.imag)))
    omega = acb(-1,arb(3).sqrt())/2
    for j in range(19):
        boundary = [target[i,j] for i in (8,7,6)]
        strength = [abs(target[14,j]+acb(-2,-5)*sum(boundary))**2/88]
        strength.extend(abs(boundary[0]+o*boundary[1]+o**2*boundary[2])**2/3
                        for o in (omega,omega.conjugate()))
        active, = [k for k,p in enumerate(strength) if p>0]
        assert all(p.contains(0) for k,p in enumerate(strength) if k!=active)
        z = max(0,int((SCALE*strength[active]).lower().floor().unique_fmpz()))
        assert arb(z)/SCALE <= strength[active]
        rows.append(("ray",active,j+1,z))
    return "".join("\t".join(map(str,row))+"\n" for row in rows)


if __name__=="__main__":
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--write",action="store_true",help="写出已认证的有限外包表")
    write = args.parse_args().write
    root = Path(__file__).resolve().parent
    path = root/FILENAME
    data = build(root)
    if write:
        path.write_text(data,encoding="utf-8")
    elif path.read_text(encoding="utf-8")!=data:
        raise ValueError("与已冻结的有限外包表不符")
    print("已写出有限外包表" if write else "全部格点与射线逐项复算一致")
