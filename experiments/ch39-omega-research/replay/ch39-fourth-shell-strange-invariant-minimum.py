"""第三步保持七维循环不变空间时，精确夹住可达的最小奇夸克质量。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def spectrum(h):
    ev = h.eig(algorithm="rump")
    assert all(v.imag.contains(0) for v in ev)
    levels = sorted((v.real for v in ev),key=lambda z:z.mid())
    assert all(x<y for x,y in zip(levels,levels[1:]))
    return levels


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    x,y,*_ = lib("fourth-shell-two-gravity-rules")["shell_data"](root)
    _,labels,_ = lib("common-embedding-central-right-ckm")["verify_nodes"](
        root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    w = labels[n+1,2]
    assert all(value >= 0 for value in w)
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    s,t = [spec(shell,2,v)[1] for v in (n+1,n+4)]
    ov = adj(t)*y*x*s
    raw = acb_mat([[ov[i,j]*w[j] for j in range(19)] for i in range(19)])*adj(ov)
    sites = shell["geo"]["sites"]
    cycle = [sites.index((b,c,a)) for a,b,c in sites]
    groups = lib("charge-six-projector-algebra-transport")["groups"](t,cycle)
    def restricted(indices):
        return spectrum(acb_mat([[raw[i,j] for j in indices] for i in indices]))
    vals,other = restricted(groups[0]),restricted(groups[1]+groups[2])
    r3 = arb(3).sqrt()
    # 两字符的射线强度在十二维补空间中各有六项为零。
    omega = acb(-arb(1)/2,r3/2)
    for phase in (omega,omega.conjugate()):
        probs = [abs(t[8,j]+phase*t[7,j]+phase**2*t[6,j])**2/3 for j in groups[1]+groups[2]]
        probs = sorted((max(arb(0),v.lower()) for v in probs),reverse=True)
        floor = sum(a*b for a,b in zip(probs,other,strict=True))
        assert floor>arb(".0004")
        print("非平凡字符的统一下界",floor.lower().str(14),flush=True)
    support = groups[0]
    rays = acb_mat([[t[14,j] for j in support],[sum(t[k,j] for k in (8,7,6))/r3 for j in support]])
    order = (2,5,3,4,6,1,0)
    chosen = acb_mat([[rays[i,j]*vals[order[j]] for j in range(7)] for i in range(2)])*adj(rays)
    ev,vectors = chosen.eig(right=True,algorithm="rump")
    k = max(range(2),key=lambda j:ev[j].real.mid())
    assert ev[k].real>ev[1-k].real
    v = acb_mat([[vectors[i,k]] for i in range(2)])
    v = v*(1/(adj(v)*v)[0,0].real.sqrt())
    strengths = [abs(z)**2 for z in (adj(v)*rays).entries()]
    seq = sorted(range(7),key=lambda j:order[j])
    assert all(strengths[i]>strengths[j] for i,j in zip(seq,seq[1:]))
    optimum = ev[k].real
    assert arb(".000276228615") < optimum < arb(".000276228617")
    second = arb(".000271")
    distance = (optimum.sqrt()-second.sqrt())/max(a.upper() for a in w).sqrt()
    leakage = distance**2*(1-distance**2/4)
    assert distance > arb(".00029159") and leakage > arb("8.50270e-8")
    print("最小奇质量、达窗所需的算子距离及泄漏概率下界",
          optimum.str(15),distance.lower().str(12),leakage.lower().str(12),flush=True)
    return optimum


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
