"""把 23226 次输入逐步放在带弱荷的费米子环上。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
from sympy import chebyshevt, expand, symbols


def verify():
    # 检查在名义相同的范数和权重下，长环的二次递推。
    L, r = symbols("L r", nonzero=True)
    for size in (3, 4, 5, 8):
        twice_chebyshev = 2*r**size*chebyshevt(size, (L+r*r/L)/(2*r))
        assert expand(twice_chebyshev - L**size - (r*r/L)**size) == 0
    ctx.prec = 150
    winding = run_path(str(Path(__file__).with_name("ch39-first-arithmetic-plaquette-winding.py")))
    size = 23226
    def alignment(n):
        return (size*winding["angle"](n)-2*arb.pi()/3).cos()
    best, second = alignment(11431), alignment(1144)
    assert best > second > 0
    factor = 2*(best-second)/(8*arb.pi()**2)
    small = factor*2*(arb(1)/2)**size/((size-1)*(size-2))
    equal_lower = factor*6/((size**2-1)*(size**2-4))
    equal_upper = 4*equal_lower
    assert small < arb("2.71e-7011")
    assert arb("7.95e-28") < equal_lower < equal_upper < arb("3.181e-27")
    beta = arb(size-16)/3
    site_bound = (8*arb.pi()**2/(beta*arb(10).log())).sqrt()
    assert arb(".06657") < site_bound < arb(".06658")
    print("23226 步带荷环的能量差 / M^4，半幅上界", small.str(12),
          "等幅区间", equal_lower.str(12), equal_upper.str(12),
          "全部开启后一数量级的站点耦合上界", site_bound.str(12))


if __name__ == "__main__":
    verify()
