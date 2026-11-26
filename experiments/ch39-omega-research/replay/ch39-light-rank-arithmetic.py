"""五种荷的轻模数作为标准逆高度平方权的非负读出。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


scale = run_path(str(Path(__file__).with_name("ch39-scale-and-zero-weights.py")))
mask = scale["local"]["mask"]
charges = (1,2,3,4,6)


def ranks_and_mean(total):
    levels = ((15-arb(189).sqrt())/2,(15+arb(189).sqrt())/2,(15+arb(261).sqrt())/2)
    inertia = {q:[[mask["window"]["supplier"](q,t+s*arb("3e-8")) for t in levels]
                  for s in (-1,1)] for q in charges}
    rows = []
    partial = {q:arb(0) for q in charges}
    mass = arb(0)
    for n,z in enumerate(acb.zeta_zeros(1,4096),1):
        phase = mask["window"]["geometry"]["phases"](z.imag)
        row = {}
        weight = 1/(total*z.imag**2)
        mass += weight
        for q in charges:
            low,high = [[f(phase) for f in funcs] for funcs in inertia[q]]
            assert low==high
            row[q] = 1+low[0]-low[1]+low[2]
            partial[q] += weight*row[q]
        if n<=5:
            rows.append([row[q] for q in charges])
    assert rows == [[3,3,3,2,2],[5,2,3,2,2],[3,2,3,2,2],[3,2,3,3,2],[3,2,5,2,2]]
    assert partial[6]>arb("2.00782") and mass<1
    means = {q:partial[q].union(partial[q]+19*(1-mass)) for q in charges}
    return partial[6],means


def verify():
    ctx.prec,ctx.cap=1536,3
    total = mask["time_data"]()[0]
    times = scale["clock"]["times"]()
    windows = scale["local"]["required_windows"](times)
    weights = scale["small_scale_weights"](times,windows[0][1],arb(2)**-176)
    standard = [1/(total*z.imag**2) for z in acb.zeta_zeros(1,5)]
    assert standard[0]/19 > arb(".0114")
    ratios = [w/p for w,p in zip(weights[:5],standard,strict=True)]
    ctx.prec=128
    mean,means = ranks_and_mean(total)
    r1,r2,r3,r4,r5 = ratios
    other = [(r2-r3)/2, r1-r3, (r5-r3)/2, r4-r3]
    other_cap = sum(v.upper() for v in other)
    d6 = ((1-r1).upper()+3*other_cap)/(mean-2)
    fraction = (1-r1).upper()+3*other_cap+2*d6
    assert fraction < arb("7.8e-9")
    print("荷六前 4096 枚对全谱平均轻模数的贡献",mean.str(13))
    print("可行模型中五荷轻模计数项占全部算术权的上界",fraction.str(12))
    # 同时纳入重模时允许明暗响应差取两种符号。
    slopes = dict(zip(charges[:4],other,strict=True))
    first_ranks = dict(zip(charges,(3,3,3,2,2),strict=True))
    residual = 1-r1-sum(slope*(means[q]-first_ranks[q]) for q,slope in slopes.items())
    slopes[6] = residual/(means[6]-2)
    total_bias = 19*sum(abs(slope) for slope in slopes.values())
    print("十九倍轻重效率差绝对值之和，按平均总响应归一",total_bias.upper().str(10))
    assert total_bias < arb("2.212e-7")


if __name__ == "__main__":
    verify()
