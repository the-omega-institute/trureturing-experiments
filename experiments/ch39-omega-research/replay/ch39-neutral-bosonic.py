"""七条件交点处的完整 A33−AQ3，一圈背景 Feynman 规范。"""
from pathlib import Path
from runpy import run_path
import sys
from flint import acb, arb, ctx
import sympy as sp


def hessian_charges():
    # 复双重态为 (phi+,(v+H+iχ)/sqrt(2))，用四个规范化实分量。
    I=sp.I
    G=[sp.Matrix([[0,1],[1,0]])/2,sp.Matrix([[0,-I],[I,0]])/2,sp.diag(1,-1)/2]
    Q=sp.diag(1,0)
    def real_gen(T):
        m=-I*T
        return m.as_real_imag()[0].row_join(-m.as_real_imag()[1]).col_join(
            m.as_real_imag()[1].row_join(m.as_real_imag()[0]))
    h=sp.Matrix([0,1,0,0])
    T=real_gen(G[2])
    assert real_gen(Q)*h==sp.zeros(4,1)
    for j,gen in enumerate((*G[:2],G[2]-sp.Symbol("s",real=True)*Q,Q)):
        coupling=2*real_gen(gen)*T*h
        support=[i for i in range(4) if coupling[i]!=0]
        assert support==([0],[2],[1],[])[j]
        assert coupling.dot(coupling)==(sp.Rational(1,4) if j<3 else 0)
    # 复带电 Goldstone 混合差权 −1/4；中性实对 +1/4。
    assert (T.T*T)[0,0]-(real_gen(Q).T*T)[0,0]==-sp.Rational(1,4)
    assert (T.T*T)[1,1]==sp.Rational(1,4)


def j_neutral(r):
    return r*((1+r)*(1/r).log1p()-1)


def scalar_difference(w,z):
    # ∫h[ℓ(z/(1-u))−ℓ(w/h)]du；解析分出端点 h log(u)。
    def integrand(u,analytic):
        f=1-u+z
        g=u*(1-u)+w
        return u*f*f.log(analytic=analytic)-g*g.log(analytic=analytic)
    value=acb.integral(integrand,0,1,abs_tol=arb("1e-33"),rel_tol=arb("1e-33"))
    assert value.is_finite() and value.imag.contains(0)
    return value.real-z*z.log()/2+w*w.log()-arb(5)/36


def verify(root,axial):
    hessian_charges()
    ctx.prec=224
    prev=run_path(str(axial))
    energy=arb("1.1763969e14").union(arb("1.17639745e14"))
    portal=arb(".42995401871").union(arb(".42995401873"))
    system=run_path(str(root/"ch39-six-broken-conditions.py"))["broken_system"](root,energy,portal)
    weights,_,det=system["weights_at"](portal)
    assert not det.contains(0) and sum(weights).contains(1)
    data=system["original"]
    rho=system["rho"]
    zw=system["weak_mass"]
    zz=(arb(".63811")**2+3*arb(".467453")**2/5)*rho/2
    zt=arb(".8616")**2*rho
    S=sum(abs(c).upper()/x for c,x in zip(weights,data["sites"],strict=True))
    p0,p1=data["p0"],data["p1"]
    def change(function,z):
        return prev["endpoint"](function,z/p0)-prev["endpoint"](function,z/p1)
    pi=arb.pi()
    # 扣除零动量后，两个块外环都由相同的 E1 增量控制。
    w_mix=-(2*change(prev["J"],zw)+arb(0,(2*S*zw/6).upper()))/(4*pi)
    z_mix=-(change(j_neutral,zz)+arb(0,(S*zz/6).upper()))/(4*pi)
    low=scalar_difference(zw/p0,zz/p0)
    # 0≤hℓ(r/h)≤r[1+log(1+1/(4r))]。高端各积分正，取其差的区间。
    def small_bound(z):
        r=(z/p1).upper()
        return (r*(1+(1/(4*r)).log1p())).upper()
    high=(-small_bound(zw)).union(small_bound(zz))
    scalar=(low-high+arb(0,(S*(zw+zz)/6).upper()))/(8*pi)
    top=3*(change(prev["J"],zt)+arb(0,(S*zt/6).upper()))/(8*pi)
    top-=(change(prev["L"],zt)+arb(0,(S*zt).upper()))/(24*pi)
    boson=w_mix+z_mix+scalar
    total=boson+top
    for name,result in (("带电 W 与 Goldstone 两次块外插入",w_mix),("中性 Z 与径向 Higgs 两次块外插入",z_mix),
                        ("中性不等质量标量对减去带电 Goldstone 标量迹",scalar),("全玻色余量",boson),("顶场余量",top),
                        ("第八条件 A33−AQ3 及七式根处 A33−A色",total)):
        print(name,result.lower().str(23),result.upper().str(23),flush=True)
    assert arb("-.074256577194663365546")<total<arb("-.074256577194663365544")
    # 从原始角核以单调上下矩形复核低端标量积分，包含两个零质量端点。
    def entropy(h,b):
        if h==0 or b==0:
            return arb(0)
        return h*(b/h).log1p()+b*(h/b).log1p()
    raw=arb(0)
    n=1024
    for k in range(n):
        a=arb(k)/n
        b=arb(k+1)/n
        left,right=a*(1-a),b*(1-b)
        lo=min(left,right)
        hi=max(left,right)
        neutral=entropy(lo,a*zz/p0).union(entropy(hi,b*zz/p0))
        charged=entropy(lo,zw/p0).union(entropy(hi,zw/p0))
        raw+=(neutral-charged)/n
    assert raw.is_finite() and raw.contains(low) and raw.rad()<arb(".01")
    print("独立的原角核分片积分与解析积分",raw.str(9),low.str(16))


if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root,root/"ch39-neutral-axial-top.py")
