"""单位门户劈裂下，非负固有时轮廓精确抵消一个弱三分量差。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    ctx.prec=224
    mp,vev,mz,energy=map(arb,("1.22089e19","246.22","91.1876","1.176397e14"))
    g1,g2,yt=map(arb,(".467453",".63811",".8616"))
    r=(vev/mp)**2/2
    w,z,t=(r*k for k in (g2**2/2,(g2**2+3*g1**2/5)/2,yt**2))
    p0,p1=(mz/mp)**2,(energy/mp)**2
    delta=p1-p0
    signs,masses=run_path(str(root/"ch39-rotation-portal-top-energy.py"))["rotation_signs"](root)
    portal=arb(".42995401872")
    centers=[a+s*portal*r for a,s in zip(masses,signs,strict=True)]
    assert len(centers)==16 and min(centers)>r
    pi=arb.pi()
    # 带电分量比原值轻 r，中性分量保留原值。所有角矩均保留有限端点。
    def scaled(tau):
        # 矩形段上 tau 为区间。积分测度 d tau / tau；返回 k/(delta*r*tau**2)。
        low=arb(0)
        high=arb(0)
        angular=arb(0)
        jcount=64
        for j in range(jcount):
            a,b=arb(j)/jcount,arb(j+1)/jcount
            sl=a.union(b)
            h=((1-a)*a).union((1-b)*b)
            if a<=arb(1)/2<=b:
                h=h.union(arb(1)/4)
            f=h*h*((-p1*h*tau).exp().union(1))/2
            fp=h*((-p1*h*tau).exp().union(1))
            ew,ez,et=(-w*tau).exp(),(-sl*z*tau).exp(),(-t*tau).exp()
            # (e^{-sz tau}-e^{-w tau})/(r tau)，均值定理避免消去。
            diff=(w-sl*z)/r*(-sl*z*tau).union(-w*tau).exp()
            block=4*w/r*ew+2*z/r*ez-3*t/r*et+diff
            top=t/r*(-t*tau).exp().union(1)
            boson_top=f*block+2*h*fp*top
            low-=boson_top/(8*pi*jcount)
            weight=3*(1-2*sl)**2
            angular+=weight*fp/jcount
        for center in centers:
            mass_delta=(-center*tau).exp().union((-(center-r)*tau).exp())
            high+=mass_delta*angular/(48*pi)
        return low+high
    # 固有时 1 至 2 给出严格正面积，且单个 Planck 节点的总残差为负。
    short=scaled(arb(1).union(2))
    area=short*delta*r*arb(3)/2
    print("短时无量纲密度、1 至 2 的弱差面积",short.lower().str(16),short.upper().str(16),
          area.lower().str(16),area.upper().str(16),flush=True)
    assert arb("0.00781")<short<arb(".01054")
    assert arb("2.2140e-46")<area<arb("2.9836e-46")

    # tau>=4/z 时，在小角段 [1/(z*tau),2/(z*tau)] 有 h*tau>=1/(2*z)。
    tail=4/z
    length_at_node=min(centers)-r
    assert length_at_node*tail>3
    time_limit=1/(2*z)
    # 将 F 的正积分截于 v<=1/p1。
    cutoff=min(time_limit,1/p1)
    fstar=delta*(-p1*cutoff).exp()*cutoff**2/2
    # 负项的简单长时下界：e^{-s z tau} >= e^{-2} 且 2z tau。
    # 同一 F 上其余轻负项只剩轴流顶；以 3t exp(-t tau) 对中性主项比较。
    assert 2*z/(3*t)*((t-z/2)*tail).exp()>2
    assert w*tail>1
    major=(16/(48*pi))*(-length_at_node*tail).exp()*tail
    # k_light <= -exp(-2)*fstar/(8*pi*tau)；重项 <=16 exp(-m tau)/(48pi)。
    assert t>z/2 and major<(-arb(2)).exp()*fstar/(8*pi)
    print("保证全尾部负号的固有时起点",tail.upper().str(15),flush=True)
    # 此处 k_light<k_total<0，所以只须用三角不等式控制原轻场的尾积分。
    # F<=1/p0，a*exp(-a)<=1/e，角积分 exp(-s z tau)<1/(z tau)。
    stop=arb("1e100")
    assert stop>tail
    tail_abs=(13/arb(1).exp()+3+1/(z*stop))/(8*pi*p0*stop)
    amplitude=tail_abs/area.lower()
    assert amplitude < arb("2.508e-21")
    print("三原子归一轮廓短时高度的严格上界",amplitude.upper().str(16),flush=True)

if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
