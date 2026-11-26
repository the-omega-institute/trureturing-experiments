"""在联合编号起的七个实际高度构造两组精确相等的双频读数。"""
from flint import acb, arb, arb_mat, ctx


N = 2_306_888_100_085


def verify():
    ctx.prec = 384
    heights = [z.imag for z in acb.zeta_zeros(N, 8)]
    times = [h-heights[0] for h in heights]
    assert all(a < b for a, b in zip(times, times[1:])) and times[-1] < 2
    reference = (arb(745)/867).acos()
    frequencies = ((arb(505)/867).acos()-reference,
                   (arb(465)/867).acos()-reference, arb(1), arb(3)/2)
    assert all(a < b for a, b in zip(frequencies, frequencies[1:]))
    assert 0 < frequencies[0] and frequencies[-1] < arb.pi()/2
    rows = [[trig for w in frequencies for trig in ((w*t).cos(), (w*t).sin())]
            for t in times]
    cofactors = []
    for j in range(8):
        minor = arb_mat([[x for k, x in enumerate(row) if k != j] for row in rows[:7]])
        cofactors.append((-1)**j*minor.det())
    assert not cofactors[0].contains(0)
    coefficients = [x/cofactors[0] for x in cofactors]
    assert all(a*a+b*b > 0 for a,b in zip(coefficients[::2],coefficients[1::2]))
    print("两组的中心化余弦、正弦系数",
          [x.str(16) for x in coefficients[:4]], [(-x).str(16) for x in coefficients[4:]])
    first = coefficients[0]+coefficients[2]
    assert arb("-.085349") < first < arb("-.085348")
    x = arb(867).sqrt()/(510**2+arb(867)).sqrt()
    y = -2*arb(867).sqrt()/(1173**2+4*arb(867)).sqrt()
    original = -36*(x*(frequencies[0]*heights[0]).sin()+y*(frequencies[1]*heights[0]).sin())
    assert arb(".27215586483699") < original < arb(".27215586483700")
    print("共同的第一项与原壳谱的第一项",first.str(19),original.str(19))
    for row in rows[:7]:
        difference = sum(a*b for a,b in zip(row,coefficients))
        assert difference.contains(0) and difference.rad() < arb("1e-65")
    gap = sum(a*b for a,b in zip(rows[7],coefficients))
    assert arb("-6.529991e-6") < gap < arb("-6.529990e-6")
    print("第八点两组读数之差",gap.str(19))


if __name__ == "__main__":
    verify()
