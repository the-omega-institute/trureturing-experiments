"""非负固有时轮廓的弱三分量差在任意高于 M_Z 的端点严格为负。"""
from flint import arb, ctx


def verify():
    ctx.prec = 160
    g1, g2, yt = map(arb, (".467453", ".63811", ".8616"))
    mz, vev = map(arb, ("91.1876", "246.22"))
    z = 1 + 3*g1**2/(5*g2**2)
    t = 2*yt**2/g2**2
    d = (mz/(g2*vev))**2
    def bound(a):
        return ((4*a-1)*(-a).exp()+(2*z*a+1)*(-z*a/2).exp()
                -3*t*a*(-t*a).exp()+2*(-(-t*a).expm1())*(-d*a).exp())
    def slope(a):
        return ((5-4*a)*(-a).exp()+(3*z/2-z*z*a)*(-z*a/2).exp()
                -3*t*(1-t*a)*(-t*a).exp()
                +2*((t+d)*(-(t+d)*a).exp()-d*(-d*a).exp()))
    # H(0)=0；先对零附近作导数下界，随后只需要函数值的下界。
    assert bound(arb(0)).is_zero()
    intervals = 256
    first = min(slope((arb(i)/intervals).union(arb(i+1)/intervals)).lower()
                for i in range(8))
    assert first > 0
    middle = min(bound((arb(i)/intervals).union(arb(i+1)/intervals)).lower()
                 for i in range(8,intervals))
    assert middle > 0
    # a>=1 时，正项 2za exp(-za/2) 大于负项 3ta exp(-ta) 的绝对值。
    assert t>z/2
    tail = 2*z/(3*t)*(t-z/2).exp()
    assert tail>1
    print("质量比 z/w、t/w 与 p0/(4w)",z.str(12),t.str(12),d.str(12))
    print("零附近导数下界、中段函数下界、长时正负项比下界",
          first.lower().str(12),middle.str(12),tail.lower().str(12))


if __name__ == "__main__":
    verify()
