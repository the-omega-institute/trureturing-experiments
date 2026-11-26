"""全部正零点的归一行列式何时选中参考夸克混合。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_series, ctx


invariant=run_path(str(Path(__file__).with_name("ch39-flavor-invariants.py")))


def scale_interval():
    return arb("100339.66931659359210804260816"),arb("100339.66931659359210804260818")


def log_xi(x):
    """返回 log ξ(1/2+x) 的一、二阶导数。"""
    s=arb_series([x+arb(1)/2,1],prec=3)
    if s[0]>12:
        # n>=2 上 log(n)^k <= n (k=1,2)，故两阶 ζ 导数绝对值均小于此尾界。
        tail=arb(2)**(4-s[0])
        z1=arb(0,tail)
        z2=arb(0,tail+tail**2)
        base=(s/2).lgamma()+s.log()+(s-1).log()
    else:
        regular=1+(s-1)*s.zeta(deflate=1)
        base=(s/2).lgamma()+s.log()+regular.log()
        z1=z2=arb(0)
    return base[1]-arb.pi().log()/2+z1,2*base[2]+z2


def verify():
    ctx.prec,ctx.cap=320,3
    c=invariant["flavor"]["singlet"](11661)[1]
    e1,e2,e3,_=invariant["invariants"]()
    masses=[arb(v) for v in (".0000135",".000268",".01401")]
    slopes=[]
    for y in masses:
        lam=y*y
        slope=(e1.derivative()(c)*lam**2-e2.derivative()(c)*lam+e3.derivative()(c))
        slope/=2*y*(3*lam**2-2*e1(c)*lam+e2(c))
        slopes.append(slope)
    assert slopes[0]>0 and slopes[1]>0>slopes[2] and sum(slopes[1:])>0
    def verify_sign(lo,hi,sign):
        # ξ 的正 theta 核使其对数导数单调，逐端点可包围整个尺度段。
        stack=[(arb(lo),arb(hi))]
        while stack:
            low,high=stack.pop()
            value=sum(d*log_xi(low*y)[0].union(log_xi(high*y)[0])
                      for d,y in zip(slopes,masses,strict=True))
            if sign*value>0:
                continue
            assert high-low>arb(".00001")
            mid=(low+high)/2
            stack.extend(((low,mid),(mid,high)))
    a,b,d=[p.derivative()(c) for p in (e1,e2,e3)]
    t=(-b-(b*b-4*a*d).sqrt())/(2*a)
    gamma=acb.zeta_zero(1).imag
    assert gamma/t.sqrt()>arb(3800)
    verify_sign(3800,100000,-1)
    verify_sign(101000,10**6,1)
    for n in range(100000,101000,5):
        scale=arb(n).union(arb(n+5))
        derivative=sum(d*y*log_xi(scale*y)[1] for d,y in zip(slopes,masses,strict=True))
        assert derivative>0

    # 大尺度仅用中、重两条谱线和 L(x) 的显式对数界。
    A=arb(10)**6
    assert A*masses[1]>100
    log_slope = slopes[1]*((A*masses[1]/(2*arb.pi())).log()/2-arb(".05"))
    log_slope += slopes[2]*((A*masses[2]/(2*arb.pi())).log()/2+arb(".05"))
    assert log_slope>0
    x=arb(100)
    assert 1/(x-arb(".5"))+1/(x+arb(".5"))+arb(".25")/x<arb(".05")
    assert 1/x+arb(2)**(4-x)<arb(".05")

    # 高精度定根时保留 ζ 项本身，避免单侧粗尾界遮住符号。
    def fine_score(scale,*,derivative=False):
        result=arb(0)
        for y,d in zip(masses,slopes,strict=True):
            s=arb_series([scale*y+arb(1)/2,1],prec=3)
            log=(s/2).lgamma()+s.log()+(s-1).log()+s.zeta().log()
            result+=d*(2*y*log[2] if derivative else log[1]-arb.pi().log()/2)
        return result
    lower,upper=scale_interval()
    assert fine_score(lower)<0<fine_score(upper)
    # 此尺度上每个零点行列式作为 c 的函数均严格凹。
    min_t=(gamma/upper)**2
    assert e2.derivative().derivative()[0]<0
    assert (e2.derivative().derivative()(arb(0,1))*min_t+
            e3.derivative().derivative()(arb(0,1)))<0
    print("完整正零点行列式的唯一共同尺度",lower.str(15),upper.str(15))
    critical=lower.union(upper)
    tail_rate=-critical*e1[1]/(2*arb.pi()*fine_score(critical,derivative=True))
    assert arb("4.2338e7")<tail_rate<arb("4.2340e7")
    for count in (1,8,64,512):
        heights=[z.imag for z in acb.zeta_zeros(1,count)]
        def finite_score(scale):
            return sum(d*y/(g*g+(scale*y)**2)
                       for d,y in zip(slopes,masses,strict=True) for g in heights)
        lo,hi=arb(3800),upper
        assert finite_score(lo)<0<finite_score(hi)
        for _ in range(36):
            mid=(lo+hi)/2
            value=finite_score(mid)
            assert not value.contains(0)
            if value>0:
                hi=mid
            else:
                lo=mid
        print("前缀枚数及其中一个驻点尺度",count,lo.union(hi).str(12))
    print("按截断高度的首阶误差系数",tail_rate.str(12))


if __name__=="__main__":
    verify()
