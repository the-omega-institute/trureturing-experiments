"""在同一个减法尺度比较两枚规范候选的共同门户。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).parent
portal = run_path(str(ROOT / "ch39-chiral-gauge-portal-mass.py"))
fixed = portal["fixed"]
chiral = fixed["chiral"]
two = fixed["two"]


def portal_at(model, initial, anchor, scale, steps=4096):
    layers, inverse, table = anchor
    unit = arb("2.435e18")/724**2
    L = (scale/model["mz"]).log()
    possible = [r for r in table if r[0]<L]
    assert all(r[0]<L or r[0]>L for r in table)
    segments = fixed["gauge_segments"](L,possible,steps)
    _,end = two["flow"]([*initial[:2],inverse],segments,fixed["sm"]["D"],fixed["sm"]["C"])
    g2,gy = 4*arb.pi()/end[1],12*arb.pi()/(5*end[0])
    subset = [unit*h for h in layers[3][3:] if unit*h<scale]
    assert len(subset)==14
    loop = 16*arb.pi()**2
    double_log = (9*g2**2+3*gy**2)*sum(m*m*(scale/m).log()**2 for m in subset)/(2*loop**2)
    correction = portal["finite_thresholds"](model,layers,unit,scale,g2,gy,double_log,report=False)
    # 80 仅是矩的数值归一化，不影响非负系数拟合的最低次数。
    moments = [sum(m*m*(1+2*(scale/m).log())*(m/(80*unit))**(2*k) for m in subset)/loop
               for k in range(5)]
    return correction,moments


def verify():
    model,_,_,initial,_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    anchors = []
    unit = arb("2.435e18")/724**2
    for n in (11661,11934):
        rows = read(acb.zeta_zero(n).imag)
        layers = {q:chiral["ordered"]([h for h,z in rows if z==q]) for q in (0,*model["charges"])}
        mu,_,alpha = fixed["crossing_enclosure"](model,initial,layers,False,724,"29.040",".1184")
        assert arb("3.72e14")<mu<arb("3.731e14")
        table = [((unit/model["mz"]).log()+l,w,B,f,s) for l,w,B,f,s
                 in two["thresholds"](model,layers,False)]
        anchors.append((layers,1/alpha,table))
    gaps = []
    residuals = []
    for i in range(8):
        lo = arb("3.72e14")+arb("1.1e12")*i/8
        hi = arb("3.72e14")+arb("1.1e12")*(i+1)/8
        (a,x),(b,y) = [portal_at(model,initial,anchor,lo.union(hi)) for anchor in anchors]
        wa,wb = x[0],y[0]
        gap = b/wb-a/wa
        assert gap>0, (i,gap)
        gaps.append(gap.lower())
        residuals.append((gap*wa*wb/(wa+wb)).lower())
    print("共用减法尺度窗中的门户差下界、最佳共同门户的质量残差下界/GeV²",
          min(gaps).str(10),min(residuals).str(10))

    # 在同一固定尺度，两个零点所需阈值比处于第零矩与第四矩之间。
    scale = arb("3.725e14")
    (a,x),(b,y) = [portal_at(model,initial,anchor,scale) for anchor in anchors]
    rr = b/a
    assert all(y[i]/x[i]<rr for i in range(4))
    assert rr<y[4]/x[4]
    fit = arb_mat([[x[0],x[4]],[y[0],y[4]]]).solve(arb_mat([[a],[b]])).entries()
    assert all(v>0 for v in fit)
    assert fit[0]+fit[1]*(scale/(80*unit))**8<arb(".08")
    print("共同非负谱多项式最低次数 4，常数项与四次项系数",*[v.str(7) for v in fit])


if __name__=="__main__":
    verify()
