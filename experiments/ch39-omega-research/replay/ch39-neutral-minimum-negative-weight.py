"""给定允许的截止能区，消去独立弱三分量差的锐负量。"""
from pathlib import Path
from runpy import run_path
from flint import arb, ctx


def phi(value):
    if value.is_zero():
        return arb(1)
    return -(-value).expm1()/value


def phi_interval(value):
    return phi(value.upper()).union(phi(value.lower()))


def integrated_gap(floor, w, z, t, p0, p1, panels):
    """对固定固有时密度消去小质量因子后，用矩形覆盖积分。"""
    R, T = z/w, t/w
    span = 1/floor-1
    total = arb(0)
    for i in range(panels):
        tau = (1+span*i/panels).union(1+span*(i+1)/panels)
        a = w*tau
        ew, et = (-a).exp(),(-T*a).exp()
        # e^{-sRa}-e^{-a} 的除差。积分指数介于 e^{-max(1,R)a} 与一之间。
        mean = (-max(R,arb(1))*a.upper()).exp().union(arb(1))
        axial = 2*T*phi_interval(T*a)
        for j in range(panels):
            lo, hi = arb(j)/panels, arb(j+1)/panels
            s = lo.union(hi)
            h = (lo*(1-lo)).union(hi*(1-hi))
            u, v = p0*tau*h, p1*tau*h
            scalar = 4*ew+2*R*(-s*R*a).exp()-3*T*et+(1-s*R)*mean
            angles = (phi_interval(u)-phi_interval(v))*scalar
            angles += axial*(-u).exp()*(-(-(p1-p0)*tau*h).expm1())
            total += h*angles
    return w*span*total/(8*arb.pi()*panels**2)


def verify(root):
    ctx.prec=192
    g1,g2,yt=map(arb,(".467453",".63811",".8616"))
    v,planck=arb("246.22"),arb("1.22089e19")
    rho=(v/planck)**2/2
    w,z,t=[rho*m for m in (g2**2/2,(g2**2+3*g1**2/5)/2,yt**2)]
    p0,p1=[(arb(m)/planck)**2 for m in ("91.1876","1.176397e14")]
    run_path(str(root/"ch39-neutral-positive-profile.py"))["verify"]()
    common=run_path(str(root/"ch39-full-neutral-free-cutoff.py"))
    high=-common["weak_column"](arb(1),p0,p1)
    low=-run_path(str(root/"ch39-neutral-mz-cutoff.py"))["low_weak_column"](w/p0,z/p0,t/p0,p0,p1)
    b_mz=low/(high-low)
    floor=arb("4e-15")
    gap=integrated_gap(floor,w,z,t,p0,p1,512)
    b_old=(high-gap)/gap
    print("最低截止为 Z 质量时的锐负量",b_mz.str(18))
    print("原支撑区间的弱差跨度、锐负量",gap.lower().str(14),gap.upper().str(14),
          b_old.lower().str(14),b_old.upper().str(14))
    assert arb("1.86547379546839") < b_mz < arb("1.86547379546841")
    assert arb("3.217e20") < b_old < arb("3.262e20")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
