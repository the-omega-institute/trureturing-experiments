"""在同一低能减法点对顶夸克移质，认证逆强耦合的位移和截止权的代价。"""
from flint import arb, ctx


def verify():
    ctx.prec = 256
    y, v, planck, mz = map(arb, (".8616", "246.22", "1.22089e19", "91.1876"))
    z = (y*v/planck)**2/2
    p0 = (mz/planck)**2
    x_min, x_max = arb("4e-15"), arb(1)
    scale = 6*arb.pi()
    a = z/p0
    b = 1+4*a
    g0 = a.log()+(3-b)*b.sqrt()*(1/b.sqrt()).atanh()+4*a
    reference = g0/scale
    assert arb(".14500786558226750640770") < reference < arb(".14500786558226750640771")

    min_p = (arb("1e13")/planck)**2
    high_end_error = 6*z/(scale*min_p)
    deficit = z/(scale*x_min)
    # 归一化使正部质量等于 (1 + 总变差)/2。单节点亏额均非负。
    error = high_end_error+13*deficit
    assert error < arb("2.62e-20")
    tolerance = arb(".01")
    necessary_tv = 2*(reference-high_end_error-tolerance)/deficit-1
    assert necessary_tv > arb("1.3485753452e20")
    print("顶质量平方比、参考色移、全窗总变差至多25的绝对误差上界",
          z.str(15), reference.str(24), error.upper().str(12))
    print("要使色移至多0.01，全窗必要总变差下界", necessary_tv.lower().str(18))

    for title, energy, limit in (
        ("全窗", "1e13", "1.399e20"),
        ("含轻弱矢量的交点盒", "1.4513765e14", "1.348805e20"),
        ("仅含轻顶的交点盒", "2.2709e14", "1.348669e20"),
    ):
        p = (arb(energy)/planck)**2
        separation = deficit*(1-(p0/4+z)/(2*x_min)-6*x_min/p)-z/(scale*x_max)
        assert separation > 0
        achievable_tv = 1+2*(reference-tolerance)/separation
        assert achievable_tv < arb(limit)
        print(title, "只约束色移及归一时，两端节点构造的总变差上界",
              achievable_tv.upper().str(17))


if __name__ == "__main__":
    verify()
