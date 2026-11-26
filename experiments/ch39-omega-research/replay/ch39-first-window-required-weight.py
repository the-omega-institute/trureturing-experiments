"""逐编号核对早期三窗候选，给第 8359 枚在共同尺度中的权重下界。"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_poly, ctx


root = Path(__file__).parent
test = run_path(str(root / "ch39-joint-flavor-time-prefix.py"))
model = run_path(str(root / "ch39-physical-shared-scale.py"))
hits = (8359, 12336, 19872, 22027)


def check_gap(block):
    return test["check_gap_block"](block)


def verify_candidates():
    jobs = [(n, min(768, right-n)) for left, right in zip(hits, hits[1:])
            for n in range(left+1, right, 768)]
    with ProcessPoolExecutor(max_workers=4) as pool:
        totals = tuple(map(sum, zip(*pool.map(check_gap, jobs), strict=True)))
    assert totals[0] == hits[-1]-hits[0]+1-len(hits) and totals[-1] == 0
    for n in hits:
        _, alpha, curvature = model["physical"]["observables"](n)
        assert arb(".117") < alpha < arb(".119") and arb(".998") < curvature < arb("1.002")
    print("逐编号三窗命中", hits, "各开间隙内复核数", totals)


def verify_weights():
    ctx.prec, ctx.cap = 640, 3
    shared = model["shared"]
    (e1, e2, _), reference, lambdas, slopes = shared["flavor_data"]()
    times = shared["clock"]["clock"]["times"]()
    first, second = [acb.zeta_zero(n).imag for n in hits[:2]]
    floor = arb("41312660.31094738")
    c_first = model["physical"]["joint"]["flavor"]["singlet"](8359)[1]
    assert e1.degree() == 1 and e2.degree() == 2 and e2[2] < 0
    print("去掉首枚后的必要尺度", (floor*second/first).str(16))
    assert floor*second/first > 57835930
    for upper, target in ((41312676, "2.63415e-9"), (50000000, "1.8171e-10"), (56000000, "1.1872e-11")):
        endpoint = (upper/second)**2
        score, (left, middle), terms = model["upper_polynomial"](endpoint, times, lambdas, slopes)
        x = (floor.union(arb(upper))/first)**2
        remainder = (x-left)**2*(x-middle)**2*(x-endpoint)
        remainder *= sum(d/(1+z*x) for d, z in zip(terms, lambdas, strict=True))
        for j in range(1, 5):
            remainder /= 1+j*j*x
        assert score < 0 and remainder.upper() > 0
        weight = -score/remainder.upper()
        assert weight > arb(target)
        print("尺度上限、强制首点权重和质量矩组合亏额", upper, weight.str(13),
              (-e2[2]*weight*(c_first-reference)**2).str(13))

    # 采用标准 1/γ² 再乘任意非负不增效率。将最后一枚以后的所有零点均计入上界。
    h, neighbor = [z.imag for z in acb.zeta_zeros(hits[-1], 2)]
    cutoff = ((h+neighbor)/2).mid()
    assert h < cutoff < neighbor and cutoff.zeta_nzeros().unique_fmpz() == hits[-1]
    tail = (cutoff/(2*arb.pi())).log()/(arb.pi()*cutoff)
    tail += (arb(7)/8-hits[-1]+arb(".28")*(cutoff.log()+arb(1)/2))/cutoff**2
    mass_upper = tail + sum(1/acb.zeta_zero(n).imag**2 for n in hits)
    minimal = 1/(first**2*mass_upper)
    assert minimal > arb(".0001903687")
    print("按高度衰减的首点概率下界", minimal.str(12),
          "质量矩组合亏额", (-e2[2]*minimal*(c_first-reference)**2).str(12))


    # 四项时间矩给最先零点的权重设上限。平方多项式比值在两根右侧递增。
    endpoint = (floor/first)**2
    _, (left, middle), _ = model["upper_polynomial"](endpoint, times, lambdas, slopes)
    moments = shared["clock"]["clock"]["cleared_moments"](times, 4)
    polynomial = (arb_poly([-left, 1])*arb_poly([-middle, 1]))**2
    factor = 1/polynomial(endpoint)
    for j in range(1, 5):
        factor *= 1+j*j*endpoint
    cap = factor*sum(polynomial[i]*v/200**i for i, v in enumerate(moments))
    assert 0 < cap < arb("2.634163e-9")
    ratio = minimal*(1-cap)/(cap*(1-minimal))
    assert ratio > 72200
    print("首点全尺度权重上界与相对逆平方谱的效率倍率下界", cap.str(14), ratio.str(12))


if __name__ == "__main__":
    verify_candidates()
    verify_weights()
