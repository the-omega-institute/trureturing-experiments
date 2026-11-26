"""七条件交点盒上的顶夸克双轴流横向质量阈值。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx
import sympy as sp


def spin_trace():
    # 欧氏 Hermitian Clifford 表示，把带质量项与动量项按准确矩阵代数分开核对。
    one=sp.eye(2)
    p=(sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),sp.diag(1,-1))
    gamma=[sp.kronecker_product(p[0],x) for x in p]+[sp.kronecker_product(p[1],one)]
    g5=gamma[0]*gamma[1]*gamma[2]*gamma[3]
    assert g5*g5==sp.eye(4)
    for i,g in enumerate(gamma):
        for j,h in enumerate(gamma):
            # 两个 -i 源顶点共同贡献负号。
            assert -sp.trace(g*g5*h*g5-g*h)==8*int(i==j)
            for k,k1 in enumerate(gamma):
                assert sp.trace(g*g5*k1*h*g5-g*k1*h)==0
                for k2 in gamma:
                    assert sp.trace(g*g5*k1*h*g5*k2-g*k1*h*k2)==0


def J(r):
    a=(1+4*r).sqrt()
    return 2*r*(a*(1/a).atanh()-1)


def L(r):
    a=1+4*r
    return r.log()+a-1+(3-a)*a.sqrt()*(1/a.sqrt()).atanh()


def endpoint(f,r):
    return f(r.lower()).union(f(r.upper()))


def verify(root):
    spin_trace()
    ctx.prec=224
    energy=arb("1.1763969e14").union(arb("1.17639745e14"))
    portal=arb(".42995401871").union(arb(".42995401873"))
    build=run_path(str(root/"ch39-six-broken-conditions.py"))["broken_system"]
    system=build(root,energy,portal)
    weights,_,det=system["weights_at"](portal)
    assert not det.contains(0) and sum(weights).contains(1)
    data=system["original"]
    z=arb(".8616")**2*system["rho"]
    r0,r1=z/data["p0"],z/data["p1"]
    size=sum(abs(c).upper()*z/x for c,x in zip(weights,data["sites"],strict=True))
    # E1(y)+log(y) 的导数在 [0,1]，每个端点的修正同号。
    I=J(r0)-endpoint(J,r1)+arb(0,(size/6).upper())
    mass=3*I/(8*arb.pi())
    V=(L(r0)-endpoint(L,r1)+arb(0,size.upper()))/(4*arb.pi())
    diagonal=V/2+mass
    relative=-V/6+mass
    print("顶质量平方与 M_Z² 之比、两处无截止双轴流矩",r0.str(16),J(r0).str(19),endpoint(J,r1).str(12))
    for name,value in (("双轴流额外质量项",mass),("顶对 A33 的完整破缺减对称位移",diagonal),
                       ("顶对 A33−A色（亦为 A33−AQ3）的破缺减对称位移",relative)):
        print(name,value.lower().str(24),value.upper().str(24))
    # 对角变量 t∈[0,1/2]，用正的分片区间积分独立检查 J 两端之差。
    nr=512
    direct=arb(0)
    for k in range(nr):
        u=arb(k)/(2*nr)
        v=arb(k+1)/(2*nr)
        left,right=u*(1-u),v*(1-v)
        def segment(r):
            return (r*(left/r).log1p()).union(r*(right/r).log1p())
        direct+=(segment(r0)-segment(r1))/nr
    direct=3*(direct+arb(0,(size/6).upper()))/(8*arb.pi())
    assert direct.is_finite() and direct.contains(mass) and direct.rad()<arb(".0001")
    assert arb(".0191957894453457348076")<mass<arb(".0191957894453457348080")
    assert arb(".1279516886320463646131")<diagonal<arb(".1279516886320463646141")
    assert arb("-.0170561769502211417945")<relative<arb("-.0170561769502211417938")
    print("分片积分的独立包围",direct.str(10))


if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
