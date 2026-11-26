"""同一 M_Z 减法处方下电磁与色逆耦合组合的背景位移。"""
from pathlib import Path
from runpy import run_path
import sys
from flint import arb, fmpq


def check_charge_kernels():
    # 在常数项与 a=(1-2u)² 的系数上逐项核验。
    gauge=(fmpq(12,11),-fmpq(3,11))
    scalar=(0,3)
    massive=(fmpq(8,7),-fmpq(3,7))
    assert all(-fmpq(22,3)*g+fmpq(1,3)*s == -7*w for g,s,w in zip(gauge,scalar,massive))
    assert fmpq(16,9)-fmpq(8,3)*fmpq(2,3)==0


def vector_log_moment(r):
    # 由 t=1-2u 对归一矢量角权的对数精确积分。
    a=1+4*r
    return r.log()+8*r/7+2*(8-a)*a.sqrt()/7*(1/a.sqrt()).atanh()


def scalar_displacement(masses,signs,du,data,weights,panels=1024):
    p0,p1=data["p0"],data["p1"]
    intervals=[]
    for i in range(panels):
        left=arb(i)/(2*panels)
        right=arb(i+1)/(2*panels)
        h=(left*(1-left)).union(right*(1-right))
        a=((1-2*right)**2).union((1-2*left)**2)
        intervals.append((h,3*a/panels))
    total=arb(0)
    for z,sign in zip(masses,signs,strict=True):
        dz=sign*du
        # 对质量平移使用精确均值积分，把完整质量路径包入导数。
        zpath=z.union(z+dz)
        assert zpath>0
        for x,c in zip(data["sites"],weights,strict=True):
            term=arb(0)
            for h,weight in intervals:
                a=zpath+p0*h
                b=zpath+p1*h
                difference=(-a/x).exp()/a-(-b/x).exp()/b
                term+=weight*difference
            total+=c*dz*term
    return total/(12*arb.pi())


def verify(root):
    check_charge_kernels()
    gauge=run_path(str(root/"ch39-all-field-signed-crossing.py"))
    signs,masses=run_path(str(root/"ch39-rotation-portal-top-energy.py"))["rotation_signs"](root)
    energy=arb("1.4513765e14").union(arb("1.4513786e14"))
    port=arb(".429954015").union(arb(".429954021"))
    data=gauge["joint_system"](root,energy,order=16,last=46)
    assert not data["matrix"].det().contains(0)
    weights=data["matrix"].solve(data["target"]).entries()
    assert data["matrix"][0,0]==1 and all(data["matrix"][0,j]==1 for j in range(5)) and data["target"][0,0]==1
    assert all(z.overlaps(a[2]) for z,a in zip(masses,data["rows"][49:65],strict=True))
    assert [j for j,s in enumerate(signs,1) if s<0]==[3,5,7,8,14,16]
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    zw=arb(".63811")**2*rho/2
    factor=7/(4*arb.pi())
    r0=zw/data["p0"]
    r1=zw/data["p1"]
    # 大端点比值虽极小，其对数项用 r1 上下端分别计算，利用 L(r) 单调。
    tail=vector_log_moment(r1.lower()).union(vector_log_moment(r1.upper()))
    error=sum(abs(c).upper()*zw/x for x,c in zip(data["sites"],weights,strict=True))
    w=factor*(tail-vector_log_moment(r0)+arb(0,error.upper()))
    s=scalar_displacement(masses,signs,port*rho,data,weights)
    total=w+s
    print("W 质量平方与 M_Z² 的比、低端矢量对数矩、高端对数矩、五权截止余量",
          r0.str(15),vector_log_moment(r0).str(17),tail.str(12),error.str(9))
    print("W 位移上下端", w.lower().str(24),w.upper().str(24))
    print("十六层重门户位移上下端", s.lower().str(18),s.upper().str(18))
    print("破缺 B 减对称 B 的上下端",total.lower().str(24),total.upper().str(24))
    assert arb("1.6692e-25")<s<arb("1.8348e-25")
    assert arb("-1.011855602472449987499")<total<arb("-1.011855602472449987493")
    # 用原减法积分独立复核质量平移的符号与归一化。
    rule=gauge["angular_rule"](14,41)
    ratio=(data["p1"]/data["p0"]).log()
    def exposure(mass,x):
        return gauge["angular_exposure"](2,mass,x,data["p0"],data["p1"],ratio,rule)
    direct=factor*sum(c*(exposure(zw,x)-exposure(arb(0),x)) for c,x in zip(weights,data["sites"],strict=True))
    assert direct.overlaps(w) and direct<0
    print("原有限动量核的直接 W 位移包围",direct.str(10))


if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
