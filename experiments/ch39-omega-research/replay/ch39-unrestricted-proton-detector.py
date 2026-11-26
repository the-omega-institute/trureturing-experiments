"""核验三枚实际正零点对任意逐面质子因子的固定检出余量。"""

from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


INDICES = (2_407_782,2_492_462,1_308_440)
COEFFICIENTS = (353198,256731,-390071)


def verify(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))
    types,_,_ = load("fixed-word-proton-faces")["signals"](root)
    geometry = load("mincut-real-triangle-tomography")
    model = load("mincut-global-winding")
    local = load("mincut-reference-ckm-action")
    shell = load("three-ray-chiral-boundary")
    flavor = load("chiral-neutral-flavor-phase")
    ctx.prec = 256
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    weights = load("height-independent-flavor-vacuum")["minimum_order_weights"](
        model,model["terms"](shell),x,y)
    rows = []
    for n in INDICES:
        h = acb.zeta_zero(n).imag
        delta = arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
        phases = shell["geo"]["phases"](h)
        row = []
        for (a,b),pair,soft in types:
            loop = geometry["unit"](pair)*phases[724]**a*phases[804]**b
            row.append((1-delta*soft)*sum(w*(loop**q).real for q,w in weights.items()))
        rows.append(row)
        print("实际零点、其高度及两只钟的主相角",n,h.str(20),
              phases[724].arg().str(15),phases[804].arg().str(15))
    denom = sum(map(abs,COEFFICIENTS))
    scores = [sum(s*r[i] for s,r in zip(COEFFICIENTS,rows))/denom for i in range(32)]
    floor = min(z.lower() for z in scores)
    assert floor>arb(".00257774")
    print("整数分离证书给全部逐面概率的检出下界",floor.str(18))
    # 每处钟角各坐标移动 1e-5 时，钟差移动至多 2e-5。亏额可在允许区间内变化。
    cap = arb(23)/25*arb(2)**-20/14**2
    robust = floor-2*arb("1e-5")*sum(q*w for q,w in weights.items())-cap
    assert robust>arb(".00155599")
    print("三处钟角各坐标移动 1e-5 后的共同检出余量",robust.str(16))
    support = (5,27,31)
    signed = [[(1 if u>0 else -1)*v for v in row] for u,row in zip(COEFFICIENTS,rows)]
    system = arb_mat([[row[i] for i in support]+[-1] for row in signed]+[[1,1,1,0]])
    primal = system.solve(arb_mat([[0],[0],[0],[1]])).entries()
    dual = system.transpose().solve(arb_mat([[0],[0],[0],[-1]])).entries()
    assert min(primal)>0 and min(dual[:3])>0 and (sum(dual[:3])-1).contains(0)
    assert (primal[-1]+dual[-1]).contains(0)
    for i in range(32):
        if i not in support:
            assert sum(u*row[i] for u,row in zip(dual,signed))>primal[-1]
    assert arb(".003215")<primal[-1]<arb(".003217")
    print("最坏三个轨道及其组概率和可达最大模",support,*[z.str(19) for z in primal])


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
