"""两个不同壳路径在公共弱站点的一圈规范门户。"""
from flint import arb, ctx
from sympy import Matrix, Rational, expand, symbols


def verify():
    gi,gj,gk,x,y = symbols("gi gj gk x y", real=True)
    u,v = Matrix([gi,-gj,0]), Matrix([gi,0,-gk])
    # x,y 为各边范数平方。每种生成元有相同的三个站点质量平方矩阵。
    mass = (x*u*u.T+y*v*v.T)/2
    trace = expand(3*(mass*mass).trace())
    assert trace.coeff(x,1).coeff(y,1)==Rational(3,2)*gi**4
    assert trace.coeff(x,2)==expand(Rational(3,4)*(gi**2+gj**2)**2)
    ctx.prec = 128
    coefficient = 9/(64*arb.pi()**2)
    initial = arb(".1")
    slope = arb(16)/(3*8*arb.pi()**2)
    unit_log = coefficient*initial**4/(1+slope*initial**2)
    quartic = 6*arb("1.26e-17")
    # x=g_i²，解 c x²/(1+a x)=quartic。
    x_at_equal = (quartic*slope+((quartic*slope)**2+4*coefficient*quartic).sqrt())/(2*coefficient)
    assert unit_log > arb("1.42e-6") and unit_log/quartic>arb("1.88e10")
    assert x_at_equal.sqrt() < arb(".000271")
    print("g_i=0.1 单位对数带的纯规范源积分",unit_log.str(12),
          "与原径向系数上界相等的站点耦合",x_at_equal.sqrt().str(12))


if __name__ == "__main__":
    verify()
