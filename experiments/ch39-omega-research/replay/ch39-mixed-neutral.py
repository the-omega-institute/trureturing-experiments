"""电弱背景上 Q 与 T3 混合横向动能相对色逆耦合的阈值位移。"""
from pathlib import Path
from runpy import run_path
import sys
from flint import arb, fmpq, fmpq_mat


def charge_algebra():
    half=fmpq(1,2)
    t3=fmpq_mat([[half,0],[0,-half]])
    charge=fmpq_mat([[1,0],[0,0]])
    higgs=fmpq_mat([[0],[1]])
    assert charge*higgs==fmpq_mat(2,1) and t3*higgs+half*higgs==fmpq_mat(2,1)
    assert charge*t3==t3*charge
    raising=fmpq_mat([[0,1],[0,0]])
    assert (charge*raising-raising*charge)==raising==(t3*raising-raising*t3)
    assert fmpq(4,3)*3*fmpq(2,3)*(half+0)/2==fmpq(2,3)
    assert fmpq(1,3)*charge[0,0]*t3[0,0]==fmpq(1,6)
    mixed=tuple(-fmpq(22,3)*g+fmpq(1,6)*s for g,s in zip((fmpq(12,11),-fmpq(3,11)),(0,3)))
    assert mixed==(-8,fmpq(5,2))


def mixed_log_moment(r):
    # 已包含 8−5(1−2u)²/2 的未归一权。
    a=1+4*r
    return 43*r.log()/6+5*(a-1)/3+(16-5*a/3)*a.sqrt()*(1/a.sqrt()).atanh()


def verify(root):
    charge_algebra()
    gauge=run_path(str(root/"ch39-all-field-signed-crossing.py"))
    portal_integral=run_path(str(root/"ch39-em-colour.py"))["scalar_displacement"]
    signs,masses=run_path(str(root/"ch39-rotation-portal-top-energy.py"))["rotation_signs"](root)
    energy=arb("1.4513765e14").union(arb("1.4513786e14"))
    port=arb(".429954015").union(arb(".429954021"))
    data=gauge["joint_system"](root,energy,order=16,last=46)
    matrix=data["matrix"]
    assert not matrix.det().contains(0) and all(matrix[0,j]==1 for j in range(5))
    assert data["target"][0,0]==1
    weights=matrix.solve(data["target"]).entries()
    assert all(z.overlaps(a[2]) for z,a in zip(masses,data["rows"][49:65],strict=True))
    assert [j for j,s in enumerate(signs,1) if s<0]==[3,5,7,8,14,16]
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    zw=arb(".63811")**2*rho/2
    r0=zw/data["p0"]
    r1=zw/data["p1"]
    tail=mixed_log_moment(r1.lower()).union(mixed_log_moment(r1.upper()))
    error=arb(43)/6*sum(abs(c).upper()*zw/x for x,c in zip(data["sites"],weights,strict=True))
    w=(tail-mixed_log_moment(r0)+arb(0,error.upper()))/(4*arb.pi())
    portal=portal_integral(masses,signs,port*rho,data,weights)/2
    total=w+portal
    print("混合角核的低端与高端对数矩、五权余量",mixed_log_moment(r0).str(20),tail.str(12),error.str(10))
    print("W 与带电 Goldstone 位移",w.lower().str(24),w.upper().str(24))
    print("重门户位移",portal.lower().str(18),portal.upper().str(18))
    print("破缺减对称的混合动能相对色位移",total.lower().str(24),total.upper().str(24))
    # 分别按原规范与标量核作同一有限动量积分，检验混合权符号。
    rule=gauge["angular_rule"](14,41)
    logratio=(data["p1"]/data["p0"]).log()
    def shift(kind,x):
        args=(x,data["p0"],data["p1"],logratio,rule)
        return gauge["angular_exposure"](kind,zw,*args)-gauge["angular_exposure"](kind,arb(0),*args)
    direct=sum(c*(arb(22)/3*shift(3,x)-shift(1,x)/6)
               for x,c in zip(data["sites"],weights,strict=True))/(4*arb.pi())
    print("原两种角核直接相加的独立包围",direct.str(12))
    assert direct.overlaps(w) and direct<0
    assert arb("8.3463e-26")<portal<arb("9.1740e-26")
    assert arb("-1.044905623704606945806")<total<arb("-1.044905623704606945799")


if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
