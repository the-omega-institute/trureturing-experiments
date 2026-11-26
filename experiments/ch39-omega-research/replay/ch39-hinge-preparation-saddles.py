"""三体制备与线性体积作用的形状、尺度平衡及全域鞍点障碍。"""
from flint import arb, ctx
import sympy as sp


def verify_elimination():
    t, r = sp.symbols("t r")
    h4, h6 = t * (4 - t), 4 * t - 1
    v = 3 * t - t**2 - 1
    area_match = r**2 * h6 - h4
    weight = 3 * r * h4 + 2 * h6
    tangency = v.diff(t) * weight - v * weight.diff(t)
    resultant = sp.Poly(sp.resultant(area_match, tangency, t), r)
    polynomial = sp.Poly(99*r**6 + 288*r**5 - 188*r**4 - 564*r**3
                         + 32*r**2 + 288*r + 44, r)
    assert resultant == -5 * polynomial
    assert polynomial.LC() % 7 and all(polynomial.eval(j) % 7 for j in range(7))


def open_plateau(n, k, area_squared):
    scaled = (n - 1) * (n - 2) * area_squared
    assert (k - 1) * (k - 2) < scaled < k * (k - 1)


def verify_examples():
    # 固定体积系数，在 65 槽的两个基态扇区内解 H=0 与形状方程。
    n, k4, k6 = 65, 31, 29
    p4, p6 = arb(k4) / n, arb(k6) / n
    lam = 6 + 3 * arb(3).sqrt() / 2
    c = 4 * arb(5).sqrt() * lam
    a, b, d = 30*p4 - c, 3*c - 120*p4 - 80*p6, 20*p6 - c
    t = (-b + (b*b - 4*a*d).sqrt()) / (2*a)
    assert arb(".74659") < t < arb(".74660")
    cos4 = (3*t - 2) / (2 * (1+t))
    cos6 = (3 - 2*t) / (2 * (1+t))
    d4 = 10 / ((1+t)**2 * (1-cos4**2).sqrt())
    d6 = -15 / ((1+t)**2 * (1-cos6**2).sqrt())
    hp = (2*a*t + b)/16
    length_squared = -(30*k4*d4 + 20*k6*d6) / (8*arb.pi()*hp)
    assert arb("1.2299") < length_squared.sqrt() < arb("1.2300")
    open_plateau(n, k4, length_squared*t*(4-t)/16)
    open_plateau(n, k6, length_squared*(4*t-1)/16)
    assert -4*length_squared*hp**2 < 0
    print("65 槽驻点的形状、长边平方、二阶行列式",
          t.str(12), length_squared.sqrt().str(12), (-4*length_squared*hp**2).str(12))

    # 平直例的角亏逐面为零，制备基态亦为包含相位的实作用的占据基态。
    n, k4, k6 = 3, 2, 2
    length_squared = 18*n*(2*arb(3).sqrt()*k6 - 3*k4)/(5*arb.pi()*(6*k4+k6))
    open_plateau(n, k4, 5*length_squared/36)
    open_plateau(n, k6, 5*length_squared/48)
    print("三槽平直驻点的体积系数与长边平方", (6*arb(5).sqrt()).str(12),
          length_squared.sqrt().str(12))


if __name__ == "__main__":
    ctx.prec = 192
    verify_elimination()
    verify_examples()
