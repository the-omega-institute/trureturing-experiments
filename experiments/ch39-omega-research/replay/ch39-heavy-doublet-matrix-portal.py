"""实际十九站谱上，十六层重双重态的厄米矩阵门户与弱三分量残差。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
import sympy as sp


def source_algebra():
    # 两个任意厄米层矩阵与质量本底都可不对易；两个中性源各块为标量。
    a = sp.Matrix([[1, sp.I],[-sp.I, 0]])
    b = sp.Matrix([[0, 1],[1, 0]])
    d = sp.diag(2,5)
    mass = sp.diag(d+a,d+a+b)
    q, t = sp.diag(1,1,0,0), sp.diag(1,1,-1,-1)/2
    assert a*b!=b*a and d*a!=a*d and d*b!=b*d
    assert q*mass==mass*q and t*mass==mass*t
    assert t*(t-q)/3==sp.diag(-1,-1,1,1)/12


def verify(root):
    ctx.prec = 192
    source_algebra()
    energy = arb("1.1763969e14").union(arb("1.17639745e14"))
    portal = arb(".42995401871").union(arb(".42995401873"))
    six = run_path(str(root/"ch39-six-broken-conditions.py"))
    scalar = run_path(str(root/"ch39-heavy-doublet-splitting.py"))
    system = six["broken_system"](root,energy,portal)
    weights = system["weights_at"](portal)[0]
    data, rho = system["original"], system["rho"]
    sites, masses = data["sites"], system["masses"]
    centers = [z+s*portal*rho for z,s in zip(masses,system["signs"],strict=True)]
    assert len(centers)==16 and sum(abs(c) for c in weights)<arb("24.465")
    assert all(a+2*rho<b-2*rho for a,b in zip(centers,centers[1:]))
    # 每段中 D=-F' >0 且 D'=-F''<0。便可求出两种范数约定的对易极值。
    panels = 768
    angular = []
    for i in range(panels):
        l, h = arb(i)/(2*panels), arb(i+1)/(2*panels)
        angular.append(((l*(1-l)).union(h*(1-h)),
                        3*((1-2*h)**2).union((1-2*l)**2)/panels))
    def curvature(z):
        total=arb(0)
        for h, angle in angular:
            a=z+data["p0"]*h
            b=z+data["p1"]*h
            def kernel(x):
                return ((-a/x).exp()*(1+a/x)/a**2
                        -(-b/x).exp()*(1+b/x)/b**2)
            total+=angle*sum(c*kernel(x) for c,x in zip(weights,sites,strict=True))
        return total

    component_max, operator_max = arb(0), arb(0)
    for center in centers:
        tube=center+arb(0,(2*rho).upper())
        d=scalar["slope_enclosure"](tube,sites,weights,data["p0"],data["p1"],panels)
        dd=curvature(tube)
        assert d>0 and dd>0,(d,dd)
        component_max+=2*rho*d/(48*arb.pi())
        operator_max+=rho*d/(48*arb.pi())
    assert arb("2.8713e-25") < component_max < arb("3.0498e-25")
    assert arb("1.4356e-25") < operator_max < arb("1.5250e-25")
    print("分量范数约定的最大正位移上下界",component_max.lower().str(17),component_max.upper().str(17))
    print("不变量耦合范数约定的最大正位移上下界",operator_max.lower().str(17),operator_max.upper().str(17))

    def uniform(u_common, extra):
        # 对 A,B 的谱范数约束分别为 u_common、extra；逐有序本征值走至多 extra*rho。
        minimum = [z-(1+u_common+extra)*rho for z in masses]
        assert min(minimum)>0
        return 25*sum((extra*rho/z).log1p() for z in minimum)/(48*arb.pi())
    bound=uniform(arb(1),arb(1))
    threshold=uniform(arb("1e22"),arb("1e22"))
    print("全窗单位不变量耦合界、1e22 不变量耦合界",bound.upper().str(17),threshold.upper().str(17))
    assert bound<arb("5.713e-24") and threshold<arb(".06527")
    run_path(str(root/"ch39-neutral-residual-all-cutoffs.py"))["verify"]()

if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
