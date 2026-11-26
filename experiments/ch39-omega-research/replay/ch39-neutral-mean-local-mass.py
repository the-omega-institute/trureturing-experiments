"""中性高模均值定心的两步磁壳质量及其两种一圈交点。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


ROOT = Path(__file__).parent


def inputs(zero=11431):
    model = run_path(str(ROOT/"ch39-native-center-parity-crossing.py"))
    ctx.prec = 224
    charges = (0,*model["charges"])
    levels = model["shell"]["shell_spectra"](charges,zero_index=zero)
    assert len(levels[0]) == 18 and all(len(levels[q])==19 for q in charges[1:])
    center = sum(levels[0])/18
    rows = [(v*abs(v-center),q) for q in charges for v in levels[q]]
    rows.sort(key=lambda item:float(item[0]))
    assert rows[0][0]>0 and rows[6][0]>rows[5][0]
    em,weak = arb("127.930"),arb(".23129")
    initial = [3*(1-weak)*em/5,weak*em,1/arb(".1180")]
    weights = {q:[6*w for w in row] for q,row in zip(charges[1:],model["weyl"])}
    weights[3] = [x+y for x,y in zip(weights[3],model["higgs"])]
    weights[0] = [arb(0),arb(-7),-arb(21)/2]
    return model,center,rows,initial,weights


def verify():
    model,center,rows,initial,weights=inputs()
    beta,pi,mz=model["beta"],arb.pi(),model["mz"]
    # t=log(mu/kappa)，L=log(mu/M_Z)。六项的尖阈值方程对 L、t、A_U 线性。
    light = rows[:6]
    assert [q for v,q in light]==[1,1,3,1,4,3]
    matrix = arb_mat([[beta[i],sum(weights[q][i] for v,q in light),2*pi] for i in range(3)])
    rhs = arb_mat([[2*pi*initial[i]+sum(weights[q][i]*v.log() for v,q in light)] for i in range(3)])
    L,t,common = matrix.solve(rhs).entries()
    mu=mz*L.exp()
    kappa=mu*(-t).exp()
    assert rows[5][0]*kappa<mu<rows[6][0]*kappa and rows[0][0]*kappa>mz
    assert arb("1.24e15")<mu<arb("1.25e15") and 37<common<38
    print("中性均值及尖阈值尺度、质量单位、共同逆耦合",center.str(16),
          mu.str(14),kappa.str(14),common.str(12))
    print("尖阈值下全体新增质量的首尾 / GeV",
          (rows[0][0]*kappa).str(12),(rows[-1][0]*kappa).str(12))

    lbox=arb(29).union(31)
    def offsets(t):
        # 利用谱核单调性在质量区间两个端点求值，避免反双曲函数的依赖膨胀。
        out = [2*pi*v for v in initial]
        for h,q in rows:
            mass=(h*(-t).exp())**2
            def value(fn):
                lower=fn(mass.upper(),arb(1))
                upper=fn(mass.lower(),arb(1))
                ref=(-2*lbox).exp()/(8*mass)
                return lower.union(upper)-arb(0).union(ref.upper())
            if q==0:
                v=value(lambda m,p:model["sheet"](m,p,True))
                for i in range(3):
                    out[i]-=weights[q][i]*v
            else:
                f=value(lambda m,p:model["sheet"](m,p,False))
                group=model["weyl"][model["charges"].index(q)]
                for i in range(3):
                    out[i]-=6*group[i]*f
                if q==3:
                    hval=value(model["scalar"])
                    for i in range(3):
                        out[i]-=model["higgs"][i]*hval
        return out

    direction = [beta[1]-beta[2],beta[2]-beta[0],beta[0]-beta[1]]
    def projected(t):
        return sum(a*b for a,b in zip(direction,offsets(t)))
    lo,hi=arb(".02"),arb(".045")
    sign=1 if projected(hi)>0 else -1
    assert sign*projected(lo)<0<sign*projected(hi)
    for _ in range(24):
        mid=(lo+hi)/2
        v=sign*projected(mid)
        if v>0:
            hi=mid
        elif v<0:
            lo=mid
        else:
            break
    assert hi-lo < arb("1e-7")
    t=lo.union(hi)
    out=offsets(t)
    predicted=(out[0]-out[1])/(beta[0]-beta[1])
    assert lbox.contains_interior(predicted)
    mu=mz*predicted.exp()
    kappa=mu*(-t).exp()
    common=(out[0]-beta[0]*predicted)/(2*pi)
    assert arb("1.18e15")<mu<arb("1.21e15") and arb("37.70")<common<arb("37.72")
    print("有限动量减法的交点与质量单位 / GeV、共同逆耦合",mu.str(13),kappa.str(13),common.str(12))
    print("有限减法下全体新增质量的首尾 / GeV",
          (rows[0][0]*kappa).str(12),(rows[-1][0]*kappa).str(12))


if __name__ == "__main__":
    verify()
