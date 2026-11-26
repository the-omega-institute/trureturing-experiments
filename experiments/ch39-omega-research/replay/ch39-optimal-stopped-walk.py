"""公平最近邻随机游走精确实现固定真空最优字质量的最小最坏步数。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


def word_masses(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))
    realized=load("realized-vacuum-word-masses")
    local=load("mincut-reference-ckm-action")
    winding=load("mincut-global-winding")
    shell=load("three-ray-chiral-boundary")
    flavor=load("chiral-neutral-flavor-phase")
    target=load("height-independent-flavor-vacuum")
    geometry=load("mincut-real-triangle-tomography")
    ctx.prec=320
    x=arb("-.01046")
    y=x-local["reference_angles"](flavor)[0]
    weights=target["minimum_order_weights"](winding,winding["terms"](shell),x,y)
    support=realized["SUPPORT"]
    assert tuple(sorted(weights))==support
    signs=realized["WINNER"]
    theta=geometry["unit"]((108469,109)).arg()
    a=arb_mat([[(q*j*theta).cos() for j in range(107)] for q in range(107)])
    rhs=arb_mat([[int(q==0),
                  signs[support.index(q)]*weights[q].sqrt() if q in support else 0,
                  int(q==106)] for q in range(107)])
    b=a.solve(rhs)
    n=[1]+[2]*106
    c=[b[j,0]-arb(n[j])/213 for j in range(107)]
    i,j=3,14
    det=b[i,1]*b[j,2]-b[j,1]*b[i,2]
    s=(c[j]*b[i,2]-c[i]*b[j,2])/det
    h=(c[i]*b[j,1]-c[j]*b[i,1])/det
    nu=[(b[k,0]+s*b[k,1]+h*b[k,2])/n[k] for k in range(107)]
    assert all(z>arb(".0031330080774") for z in nu)
    assert (sum(w*z for w,z in zip(n,nu))-1).contains(0)
    return nu


def verify(root):
    nu=word_masses(root)
    initial=[2*sum((k-j)*nu[k] for k in range(j+1,107)) for j in range(106)]+[arb(0)]
    assert all(v>0 for v in initial[:-1])
    d=initial
    expectation=2*sum(k*k*nu[k] for k in range(1,107))
    print("目标位置方差与任意有界实现的平均步数",expectation.str(22))
    # 独立正向传播尚未作停止选择的到达质量，并把已停路径永久保留在原位。
    active=[arb(1)]+[arb(0)]*106
    held=[arb(0)]*107
    limit=4456
    for step in range(1,limit+1):
        nxt=[arb(0)]*107
        for j in range(106):
            if d[j].is_zero():
                continue
            z=(d[abs(j-1)]+d[j+1])/2-nu[j]
            if z>0:
                nxt[j]=z
            else:
                assert z<0, (step,j,z)
        released=[arb(0)]*107
        for j in range(107):
            if d[j].is_zero():
                held[j]+=active[j]
            elif nxt[j]>0:
                assert held[j].is_zero()
                released[j]=active[j]
            else:
                assert held[j].is_zero() and active[j]>d[j]
                released[j]=d[j]
                held[j]=active[j]-d[j]
        assert released[106].is_zero()
        active=[released[1]]+[(released[j-1]+released[j+1])/2 for j in range(1,106)]+[released[105]/2]
        if all(z.is_zero() for z in nxt):
            last=[(j,z.str(24)) for j,z in enumerate(d) if z>0]
            assert step==4456 and [j for j,z in enumerate(d) if z>0]==[1]
            assert arb(".00219693514864283362126")<d[1]<arb(".00219693514864283362127")
            error=[a+b-p for a,b,p in zip(held,active,nu)]
            assert all(e.contains(0) and abs(e)<arb("1e-35") for e in error)
            assert (2*d[1]-active[0]-2*sum(active[1:])).contains(0)
            print("最后一个非零势差的时刻、位置绝对值及势差",step-1,last)
            print("贪心策略需走满全部步数的概率",(2*d[1]).str(22))
            print("正向公平传播的终端逐字质量与目标之差的绝对上界",
                  max(abs(e).upper() for e in error).str(8))
            print("首次处处为精确零的步数",step,flush=True)
            return
        d=nxt
    raise AssertionError("超过预设步数")


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
