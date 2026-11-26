"""两枚三代规范交点上，电弱规范源产生 Higgs 二次项的双对数系数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpq_mpoly_ctx

ROOT = Path(__file__).parent
fixed = run_path(str(ROOT / "ch39-chiral-clock-unit-readout.py"))
chiral = fixed["chiral"]


def gauge_trace():
    # 两个同超荷双重态，x=H†H、z=S†S、u=|H†S|²，|总 Pauli 向量|²=(x-z)²+4u。
    ring = fmpq_mpoly_ctx.get(("a","b","x","z","u","v","w"))
    a,b,x,z,u,v,w = ring.gens()
    # a,b 分别为 g2²、gY²。规范质量平方矩阵的三维块为 a(x+z)I/2。
    tr4 = (3*a**2+b**2)*(x+z)**2/4 + a*b*((x-z)**2+4*u)/2
    # 扣除两个各自背景，把 3/2 Tr M_V^4 转为 16π² beta_V。
    mixed = (tr4-(3*a**2+b**2+2*a*b)*(x**2+z**2)/4)*3/2
    lam3 = (3*a**2+b**2-2*a*b)*3/4
    lam4 = 3*a*b
    assert mixed == lam3*x*z+lam4*u
    assert 2*lam3+lam4 == (9*a**2+3*b**2)/2
    # v=M²、w=λ3、b 在此恒等式中代替 λ4。标量 Hessian 的两对简并根。
    trace_s = 2*(v+w*x)**2+2*(v+(w+b)*x)**2
    assert (trace_s-4*v**2-2*x**2*(w**2+(w+b)**2))/2 == 2*v*x*(2*w+b)


def finite_thresholds(model, layers, unit, mu, g2, gy, double_log, report=True):
    # Farina, Pappadopulo, Strumia，arXiv:1303.7244v3 式 (3)、(21)、(22)。
    # 系数改记为 V 中的 +m_H² H†H，即论文 δm² 的 −1/2。
    loop2 = (16*arb.pi()**2)**2
    scalar = fermion = scalar_ll = arb(0)
    counts = []
    floor = arb(0)
    for q in model["charges"]:
        weak = 2 if q in (1,3) else 1
        color = 3 if q in (1,2,4) else 1
        group = arb(weak*weak-1)*g2**2/4+(arb(q)/6)**2*gy**2
        factor = color*weak*group/(2*loop2)
        active = 0
        for h in layers[q][3:]:
            mass = unit*h*(32 if q%2==0 else 1)
            if mass > mu:
                continue
            assert mass < mu
            active += 1
            t = (mu/mass).log()
            term = mass**2*factor
            # c=2，每个有质量的费米子是一对 Dirac 左右分量。
            fermion += 2*term*(1+12*t)
            floor += 2*term
            if q==3:
                scalar += term*(6*t**2-4*t+arb(7)/2)
                scalar_ll += 6*term*t**2
                floor += arb(17)/6*term
        counts.append((q,active))
    assert all(32*unit*h>mu for h in layers[0])
    assert scalar_ll.overlaps(double_log)
    total = scalar+fermion
    assert scalar>0 and fermion>0 and total>floor
    ref = arb("125.25")**2/2
    if report:
        print("分荷激活 Dirac",counts,"有限标量、有限费米、总额/GeV²、同尺度下界、相对电弱参数",
              *[(v.lower().str(9),v.upper().str(9)) for v in (scalar,fermion,total,floor,total/ref)])
    return total

def verify():
    gauge_trace()
    model,_,_,initial,_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    portal_intervals = []
    sensitivities = []
    for n in (11661,11934):
        rows = read(acb.zeta_zero(n).imag)
        layers = {q:chiral["ordered"]([h for h,z in rows if z==q]) for q in (0,*model["charges"])}
        assert all(layers[q][2]<9<layers[q][3] for q in model["charges"])
        mu,common,_ = fixed["crossing_enclosure"](model,initial,layers,False,724,"29.040",".1184")
        unit = arb("2.435e18")/724**2
        masses = [unit*h for h in layers[3][3:]]
        active = [m for m in masses if m<mu]
        assert all(m<mu or mu<m for m in masses)
        assert len(active)==14
        g2 = 4*arb.pi()/common
        gy = 3*g2/5
        source = (9*g2**2+3*gy**2)/2
        loop = 16*arb.pi()**2
        curv = 2*source*sum(m*m for m in active)/loop**2
        generated = source*sum(m**2*(mu/m).log()**2 for m in active)/loop**2
        ratio = generated/(arb("125.25")**2/2)
        print(n,"共同逆耦合、门户源、二次项曲率、逐层双对数贡献/GeV²、与电弱参数之比",
              *[(v.lower().str(9),v.upper().str(9)) for v in (common,source,curv,generated,ratio)])
        assert source>arb(".48") and ratio>arb("1e19")
        total = finite_thresholds(model,layers,unit,mu,g2,gy,generated)
        assert total>arb("1e25")
        sensitivity = sum(m*m*(1+2*(mu/m).log()) for m in active)/loop
        # 单独采用 V_mix=(P/2)|H|²|S|²，所有十四个活动双重态共用 P。
        target = arb("125.25")**2/2
        portal = total/sensitivity
        width = target/sensitivity
        portal_intervals.append(portal+arb(0,width.upper()))
        sensitivities.append(sensitivity)
        print("抵消门户与电弱半宽",
              *[(v.lower().str(10),v.upper().str(10)) for v in (portal,width)])
    gap = portal_intervals[0]-portal_intervals[1]
    assert not gap.contains(0), gap
    residual = abs(gap).lower()*min(v.lower() for v in sensitivities)/2
    print("同一门户用于两枚时，至少一枚的质量平方残差下界 / GeV²",residual.str(11))


if __name__=="__main__":
    verify()
