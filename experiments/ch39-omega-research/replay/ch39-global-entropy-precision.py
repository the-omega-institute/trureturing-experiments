"""全谱与截断谱熵解的乘子、径向精度和尾误差的有限证书。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

root = Path(__file__).parent
entropy = run_path(str(root / "ch39-neutral-entropy-weights.py"))
joint = entropy["joint"]
scale = 2_903_600_000


def verify():
    rows, times, _, _, reference, cp0 = joint["columns_at"](scale, joint["eight_indices"])
    ctx.prec = 640
    basis = arb_mat(list(zip(*(r[0] for r in rows), strict=True)))
    weight = joint["solve"](rows, times)
    zmax = entropy["inverse_square_tail"](8000)
    logs = [-(r[1]/(scale**2*zmax)).log() for r in rows]
    expected = sum(v*b for v,b in zip(weight,logs,strict=True))
    deviation = [max(b,expected/v) for v,b in zip(weight,logs,strict=True)]
    inverse = basis.transpose().inv()
    norm = sum(abs(inverse[i,j])*deviation[j] for i in range(1,8) for j in range(8))
    assert norm < arb("2.729e17")

    _, _, poles, slopes = joint["model"]["shared"]["flavor_data"]()
    bound_g = sum(abs(v)/z for v,z in zip(slopes,poles,strict=True))
    bound_g_prime = sum(abs(v) for v in slopes)
    strength = arb(23)/25*arb(2)**-20
    offset = joint["replay"]["neutral"]["shell"]["soft"]["point"]["study"]["offset"]
    first = acb.zeta_zero(8359).imag
    assert offset*strength/(first**2+offset) < strength/2
    x_lipschitz = norm*max(arb(756),756*bound_g+19*bound_g_prime)
    delta_lipschitz = x_lipschitz*4*scale**2/strength
    assert delta_lipschitz < arb(2)**153
    data = joint["cp"]["two"]["spectral_data"]()
    for j in range(128):
        c = arb(-1+arb(j)/64).union(-1+arb(j+1)/64)
        up, _ = joint["cp"]["two"]["probabilities"](c, data)
        assert 0 < up[0] < arb(".002") and abs(up[1]) < arb(".001")
    c_lipschitz = norm*19*arb(".001")
    assert c_lipschitz < arb(2)**53
    print("全谱乘子一范数、辐边和六面坐标的对数权重 Lipschitz 上界",
          norm.str(12), delta_lipschitz.str(12), c_lipschitz.str(12))

    x = arb(0).union((arb(scale)/2**80)**2)
    # 中性非零谱的隙大于一，η≤1+18x。
    eta = 1+arb(0).union(18*x.upper())
    g = joint["model"]["gradient"](x,poles,slopes)
    bounds = []
    for u in (arb(0), arb(".002")):
        readings = [1,*[1/(1+j*j*x) for j in range(1,6)], eta*g, eta*(u-reference[6])]
        coefficients = basis.solve(entropy["column"](readings)).entries()
        least = min((v/p).lower() for v,p in zip(coefficients,weight,strict=True))
        bound = sum((v-least*p)*b for v,p,b in zip(coefficients,weight,logs,strict=True))
        bounds.append(bound)
    assert max(bounds) < 66965
    # H=2^(100000+4b)，log H+1≤2^(20+b)，b≥0。
    assert 100000*arb(2).log()+1 < 2**20 and 4*arb(2).log() < 1+100000*arb(2).log()
    partition_bits = (arb(66965)+arb("44.074")-arb("2.95816e-8").log())/arb(2).log()
    assert partition_bits < 97000
    assert cp0 > arb(".00003") and 38/cp0 < 2**21
    print("极高处的对偶分数上界与尾积分前因子的二进位数",
          max(bounds).str(12), partition_bits.str(12))


if __name__ == "__main__":
    verify()
