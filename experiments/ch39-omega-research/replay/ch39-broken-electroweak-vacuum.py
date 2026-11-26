"""在电弱破缺真空重算两个真实零点的热核、最少补场数与标量势。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


graded = run_path(str(Path(__file__).with_name("ch39-chiral-graded-vacuum.py")))
heat = graded["heat"]


def light_spectrum(planck):
    # (物种，净偏振，牛顿权的六倍，参考质量/GeV)。中微子先取三个无质量 Weyl。
    data = [("γ,g,ν", 12, -33, "0"), ("W", 6, -6, "80.38"),
            ("Z", 3, -3, "91.1876"), ("h", 1, 1, "125.25")]
    data.extend((name, -12, 6, mass) for name, mass in
                (("u", ".00216"), ("d", ".00467"), ("s", ".093"), ("c", "1.27"),
                 ("b", "4.18"), ("t", "172.5")))
    data.extend((name, -4, 2, mass) for name, mass in
                (("e", ".000511"), ("μ", ".10566"), ("τ", "1.77686")))
    assert sum(g for _, g, _, _ in data) == -62
    assert sum(a for _, _, a, _ in data) == 1
    return [(n, g, arb(a)/6, (arb(m)/planck)**2) for n, g, a, m in data]


def minimal_partners(g, target):
    # 新增 s 个标量、d 个 Dirac 后 g+s−4d=target。
    d = max(0, (g - target + 3) // 4)
    return target - g + 4*d, d


def light_completion(light, a):
    levels = [(g, alpha, z, -4 if z.is_zero() else 0) for _, g, alpha, z in light]
    levels += [(0, arb(0), j*a, weight) for j, weight in ((1,12), (2,-12), (3,4))]
    assert all(not z.overlaps(w) for i, (_, _, z, _) in enumerate(levels)
               for _, _, w, _ in levels[:i])
    ns = nd = 0
    rows = []
    for g, alpha, z, target in levels:
        s, d = minimal_partners(g, target)
        assert g + s - 4*d == target and (s < 4 or d == 0)
        ns += s
        nd += d
        rows.append((alpha + arb(s)/6 + arb(d)/3, z))
    assert (ns, nd) == (106, 11)
    return rows


def match(rows):
    # 破缺修复后的牛顿权在零质量层可以为负，单调性仅用于标定窗口。
    lo, hi = arb(".029"), arb(".034")
    k = graded["newton"]
    assert k(rows,lo)<1<k(rows,hi)
    lower = sum(w*(-z/(lo if w>0 else hi)).exp() for w,z in rows)
    assert lower > 190
    for _ in range(48):
        mid = (lo+hi)/2
        val = k(rows,mid)
        if val < 1:
            lo = mid
        elif val > 1:
            hi = mid
        else:
            break
    return lo.union(hi)


def copy_bound_and_low_branch(base, a, magnitude, target, planck, i0):
    # base 含一份四层谱，提出单份的独立贡献后，破缺基底的无质量牛顿权为 −9/2。
    ns = [(arb(1)/3, arb(0)), (arb(2), a), (arb(1), 2*a), (arb(2)/3, 3*a)]
    pure = base + [(-w,z) for w,z in ns]
    massless = sum(w for w,z in pure if z.is_zero())
    assert (massless+arb(9)/2).contains(0)
    # 以下粗界容许在相同质量处有数行。额外的配平四标量加一 Dirac 只增牛顿项。
    neg = sum(w for w,z in pure if w < 0)
    assert neg > -100
    min_count = magnitude/(a**2*i0/arb.pi())
    assert min_count > arb("1.344e55")
    # 再以 r=v 处的一阶、二阶导数匹配 Higgs 时，每份真空贡献至多放大三倍。
    cap_refit = (arb("7.24e-9")/planck)**2
    refit_coeff = sum(w*heat(2,z/cap_refit) for w,z in ns)
    assert cap_refit*(neg + min_count*refit_coeff/3)/(2*arb.pi()) > 1
    cap = (arb("4.18e-9")/planck)**2
    coeff = sum(w*heat(2,z/cap) for w,z in ns)
    assert cap*(neg + min_count*coeff)/(2*arb.pi()) > 1

    # 低截止下按无质量近似作整数取整，再包住整个有限谱的根。
    A = 18*arb.pi()/(magnitude+target)
    real_count = (A+27+(A*(A+54)).sqrt())/2
    count = int((real_count+arb(1)/2).floor().unique_fmpz())
    x0 = 6*arb.pi()/(arb(count)-arb(27)/2)
    band = x0*(1+arb(0,"1e-85"))
    zmin = min(z.lower() for _,z in base if z>0)
    assert zmin > 0 and zmin/band.upper() > 10**6
    # E₂(t)≤exp(−t), E₃(t)≤exp(−t)/2。对 K=1 方程与真空分别估计。
    tail = (-zmin/band.upper()).exp()
    abs_newton = sum(abs(w) for w,z in pure) + 4*count
    assert abs_newton*band.upper()*tail < arb("1e-100")
    assert arb(count)/3 + neg > 0
    delta = arb("1e-85")
    # 两端的无质量 K 恰为 1±delta，总 K 的误差严格小于 delta。
    xx = (x0*(1-delta)).union(x0*(1+delta))
    vacuum0 = arb(count)*xx**2/(2*arb.pi())
    vacuum_tail = arb(count)*xx.upper()**2*7*tail/(2*arb.pi())
    ratio = (vacuum0 + arb(0,vacuum_tail) - magnitude)/target
    assert arb("6.6e-15") < ratio-1 < arb("6.7e-15")
    print("必要份数下界、统一取整份数、低截止/eV 与目标相对偏差",
          min_count.lower().str(10), count, (xx.sqrt()*planck*10**9).str(12), (ratio-1).str(9))


def verify():
    chiral = graded["minimal"]["chiral"]
    model, *_ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 768
    source = chiral["source"]["spectral_supplier"](model)
    planck, v, mh = arb("1.22089e19"), arb("246.22"), arb("125.25")
    target = arb("1e-122")
    light = light_spectrum(planck)
    i0 = (9*arb(3).log()-12*arb(2).log())/2
    tree = -arb.pi()*mh**2*v**2/planck**4
    rho_star = target*planck**4/(8*arb.pi())
    print("希格斯树级 Δ(ΛG)、目标能量密度 GeV⁴ 与相对抵消精度",
          tree.str(10), rho_star.str(10), abs(target/tree).str(10))
    for number in (22027,28325):
        src = source(acb.zeta_zero(number).imag)
        bands = {q: chiral["ordered"]([h for h,charge in src if charge==q])
                 for q in (0,*model["charges"])}
        heavy = graded["minimal"]["fields"](bands)
        paired = graded["paired_spectrum"](heavy)
        assert paired[0][1].is_zero()
        M = 32*arb("2.435e18")/724**2*bands[2][15]
        a = (v**2/(M*planck))**2
        # 保留原隐谱时，仅用破缺谱替换原带荷和中性轻场的可见合计。
        frozen = paired + graded["low_sector"](4,a)
        frozen += [(alpha,z) for _,_,alpha,z in light] + [(-heavy[0][1],arb(0))]
        x = match(frozen)
        # E₃(z/x)−1/2 的绝对值对 x 单调递减。用端点计算，避免区间相减的消去。
        def change(z):
            if z.is_zero():
                return arb(0)
            return (heat(3,z/x.lower())-arb(1)/2).union(heat(3,z/x.upper())-arb(1)/2)
        broken = -x**2*sum(g*change(z) for _,g,_,z in light)/(4*arb.pi())
        cube = a**2*i0/arb.pi()*(1-arb(0).union(a/(x*i0)))
        assert broken+cube+tree < -arb("4e-36")
        # 从头按破缺谱配平时，上游 896 个重标量和 162 个重 Dirac 继续按层补足。
        repaired = paired[1:] + light_completion(light,a)
        x_new = match(repaired)
        copy_bound_and_low_branch(repaired,a,-tree,target,planck,i0)
        vacuum = a**2*i0/arb.pi()*(1-arb(0).union(a/(x_new*i0)))
        assert abs(vacuum/target-1)<arb(".0006")
        assert -arb("1.4e-67") < vacuum+tree < -arb("1.3e-67")
        # 只保留四层项的背景依赖，重标定二次项和四次项，使 v 与 Higgs 曲率质量固定。
        y = a / x_new
        assert y < arb("1e-58")
        # f(y)≤1, exp(−y)(1−exp(−y))²/y≤y，余项按真空共同单位计。
        unit = a**2 / arb.pi()
        residual = unit * (3*i0 - arb(0).union(10*y))
        assert (vacuum+tree < 0 and residual+tree < 0)
        assert arb("2.998") < residual/target < arb("3.001")
        # 核对 r⁸ I(r⁴/(M²Q²)) 的二阶 Euler 导数余项与标量质量上界。
        upper_mass = 4*(vacuum*planck**4/(8*arb.pi())).sqrt()/v
        assert upper_mass*10**9 < arb("4.84e-17")
        print(number, "旧场表、重配后 Q/M_P，旧隐谱电弱残差、重配后一圈比值、径向质量正曲率上界/eV",
              x.sqrt().str(10), x_new.sqrt().str(10), broken.str(10),
              (vacuum/target).str(11), (upper_mass*10**9).str(10))


if __name__ == "__main__":
    verify()
