"""单步谱读数下，几类候选稳定性泛函在双钟环面上对共同秩 {1,3,7,15} 的优胜区域。

候选：共同与单荷谱隙宽度、谱隙比；phi 在谱区间上为正的 sum phi、Z、F、S；
费米能级狄拉克海 E_mu 配由边表算出的全部站点度、站点方差及其统计量。

记号。对荷 q 与环面点 (U,V)=(u^q,v^q)，十九维磁拉普拉斯按三循环特征标
分为 7+6+6 三块；块 k 的第 i 个本征值记作标签 (k,i)。全谱的共同分界 1,3,7,15
把标签分成五个区：Z1={(0,1)}，Z23={(1,1),(2,1)}，Z47={(0,2),(0,3),(1,2),(2,2)}，
Z815={(0,4),(0,5),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)}，Z16=其余四个。

证书。
1. 逐标签带界：块 k 中 t 不是外点块 F 的本征值时，Haynsworth 惯性给出
   #{本征值<t} = #{eig F<t} + [s(t;U,V)<0]，s 为中心一维舒尔补。s 关于 V 精确
   取极，关于 U 用 Lipschitz 圆扫描，于是 s 在整张环面上的统一符号可由 Arb 判定。
2. 区平均的紧下界：块 k 内任意一组标签之和随 (X,Y) 的变化不超过
   (‖b_1‖+‖b_2‖) max(|dX|,|dY|)，b_1、b_2 为两钟在中心列上的系数向量（扰动是秩二的 [[0,a],[a*,0]]，Lidskii）。
   在自适应网格的格心用 Arb 求本征值。
3. 软辐边：L_delta = L_0 - delta K，K 为三条辐边的星图拉普拉斯，0<=K<=4，
   故每个标签只下移不超过 4 delta < 4e-8，以下统一留 1e-7 的余量。
"""
from pathlib import Path
from runpy import run_path
import time

from flint import acb, acb_mat, arb, ctx

HERE = Path(__file__).with_name
shell = run_path(str(HERE("ch39-three-ray-chiral-boundary.py")))
ranks = run_path(str(HERE("ch39-all-common-shell-ranks.py")))
CHARGES = (1, 2, 3, 4, 6)
adj = lambda a: a.conjugate().transpose()
SLACK = arb("1e-7")

ZONES = {
    "Z1": [(0, 1)],
    "Z23": [(1, 1), (2, 1)],
    "Z47": [(0, 2), (0, 3), (1, 2), (2, 2)],
    "Z815": [(0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)],
    "Z16": [(0, 6), (0, 7), (1, 6), (2, 6)],
}
ZONE_WINDOWS = {"Z1": ("0", "3.1"), "Z23": ("7.1", "8.69"), "Z47": ("9.6", "11.9"),
                "Z815": ("14.85", "18.33"), "Z16": ("18.36", "23")}

# (下界 t_l, 上界 t_u)：标签在整张环面上严格落在 (t_l, t_u) 内（delta=0）。
BANDS = {
 1: [[("1.1534","2.9938"),("10.3152","10.9997"),("10.9991","11.8918"),("15.0010","15.0022"),("16.6907","17.3025"),("19.0508","19.0899"),("22.1715","22.9560")],
     [("7.4575","8.2922"),("10.2432","11.0003"),("15.0022","15.0033"),("15.4242","15.7325"),("16.9713","18.1231"),("19.1592","19.2430")],
     [("7.3699","8.0506"),("10.1596","11.0003"),("14.9888","14.9899"),("15.4336","15.6988"),("17.0320","18.1983"),("19.0706","19.1373")]],
 2: [[("0.0412","2.0405"),("9.6088","10.9970"),("10.9967","11.7548"),("15.7558","16.2208"),("17.0012","17.0047"),("18.3706","19.0078"),("19.0024","20.0739")],
     [("7.6810","8.4492"),("10.4289","10.9998"),("15.2602","15.3270"),("17.0189","17.0230"),("17.3764","17.5658"),("19.0720","20.0950")],
     [("7.4785","8.0008"),("10.2724","10.9997"),("15.1065","15.1347"),("16.9563","16.9583"),("17.3466","17.5185"),("18.9487","19.9413")]],
 3: [[("1.2278","3.0472"),("10.3165","10.9937"),("10.9960","11.8782"),("15.0126","15.0154"),("16.6802","17.2660"),("19.0581","19.0978"),("22.1804","22.9245")],
     [("7.4983","8.4948"),("10.3194","10.9990"),("14.9952","14.9963"),("15.4559","15.7874"),("16.9272","18.1444"),("19.2715","19.3649")],
     [("7.2512","7.7762"),("10.0694","10.9982"),("14.9391","14.9415"),("15.4724","15.6781"),("17.1122","18.3001"),("19.0060","19.0517")]],
 4: [[("0.1655","2.1314"),("9.6228","10.9863"),("10.9878","11.7258"),("15.7468","16.1765"),("17.0061","17.0155"),("18.3926",None),("19.0112","20.0353")],
     [("7.7084","8.5988"),("10.5018","10.9981"),("15.3179","15.4101"),("17.0059","17.0138"),("17.5146","17.6690"),("19.1746","20.2922")],
     [("7.3241","7.7159"),("10.1940","10.9971"),("15.0096","15.0238"),("16.9125","16.9168"),("17.3985","17.5357"),("18.9155","19.9150")]],
 6: [[("0.3697","2.2800"),("9.6442","10.9677"),("10.9729","11.6768"),("15.7329","16.1090"),("17.0136","17.0289"),("18.4298",None),("19.0258","19.9748")],
     [("7.6807","8.6847"),("10.5816","10.9940"),("15.3643","15.4871"),("16.9791","16.9845"),("17.7236","17.8062"),("19.3181","20.4858")],
     [("7.1531","7.4106"),("10.1328","10.9920"),("14.8994","14.9038"),("16.8791","16.8850"),("17.4544","17.5453"),("18.8993","19.9041")]],
}


def bases():
    return list(ranks["character_bases"](shell))


def secular_blocks(q, basis_list):
    base, clock = ranks["build_family"](shell, q)
    out = []
    for basis in basis_list:
        block = adj(basis)*base*basis
        terms = [adj(basis)*t*basis for t in clock]
        center, = [j for j in range(basis.ncols()) if not basis[6, j].is_zero()]
        outer = [j for j in range(basis.ncols()) if j != center]
        # The clock terms live only in the center column.
        for t in terms:
            assert all(t[i, j].is_zero() for i in range(t.nrows()) for j in range(t.ncols()) if j != center)
        fixed = acb_mat([[block[i, j] for j in outer] for i in outer])
        levels = [v.real for v in fixed.eig(multiple=True)]
        links = [acb_mat([[t[i, center]] for i in outer]) for t in [block, *terms]]
        norms = [sum((abs(t[i, center])**2 for i in outer), arb(0)).sqrt() for t in terms]
        out.append(dict(c=block[center, center].real, F=fixed, levels=levels, links=links,
                        block=block, terms=terms, lip=norms[0]+norms[1]))
    return out


def uniform_count(blk, t, count=1024):
    """返回 (#eig F<t, sign)，sign=+1/-1 表示 s 在整张环面上统一同号，0 表示未判定。"""
    assert all(v < t or v > t for v in blk["levels"])
    n = sum(1 for v in blk["levels"] if v < t)
    shifted = acb_mat(blk["F"])
    for i in range(shifted.nrows()):
        shifted[i, i] -= t
    inv = shifted.inv()
    g = [[(adj(a)*inv*b)[0, 0] for b in blk["links"]] for a in blk["links"]]
    s = (blk["c"]-t-sum(g[i][i] for i in range(3))).real
    off = 2*sum(abs(g[i][j]) for i, j in ((0, 1), (0, 2), (1, 2)))
    if s > off:
        return n, 1
    if s < -off:
        return n, -1
    sign = 1 if s > 0 else -1
    rad = arb.pi()/count
    for i in range(count):
        u = acb(0, 2*arb.pi()*i/count).exp()
        const = s-2*(u*g[0][1]).real
        amp = 2*abs(g[0][2]+u.conjugate()*g[1][2])
        var = 2*(abs(g[0][1])+abs(g[1][2]))*rad
        if not sign*const-amp > var:
            return n, 0
    return n, sign


def certify_bands(blocks_by_q):
    lo, hi = {}, {}
    for q in CHARGES:
        for k in range(3):
            for i, (tl, tu) in enumerate(BANDS[q][k], start=1):
                n, st = uniform_count(blocks_by_q[q][k], arb(tl))
                assert n+(st <= 0) <= i-1, (q, k, i, tl)
                lo[q, k, i] = arb(tl)-SLACK
                if tu is not None:
                    n, st = uniform_count(blocks_by_q[q][k], arb(tu))
                    assert n+(st < 0) >= i, (q, k, i, tu)
                    hi[q, k, i] = arb(tu)
    # Every label sits inside its zone window, so zones are position sets.
    for q in CHARGES:
        for name, labels in ZONES.items():
            a, b = map(arb, ZONE_WINDOWS[name])
            for lab in labels:
                assert lo[(q,)+lab] > a
                if (q,)+lab in hi:
                    assert hi[(q,)+lab] < b
        assert sum(len(v) for v in ZONES.values()) == 19
    return lo, hi


def block_eigs(blocks, X, Y):
    U, V = acb(0, X).exp(), acb(0, Y).exp()
    result = []
    for blk in blocks:
        P = U*blk["terms"][0]+V*blk["terms"][1]
        M = blk["block"]+P+adj(P)
        result.append(sorted((x.real for x in M.eig(multiple=True)), key=lambda a: a.mid()))
    return result


def zone_mean(eigs, name):
    labels = ZONES[name]
    for zname, zlabels in ZONES.items():
        a, b = map(arb, ZONE_WINDOWS[zname])
        for k, i in zlabels:
            assert a < eigs[k][i-1] < b
    return sum(eigs[k][i-1] for k, i in labels)/len(labels)


def adaptive_inf_a47(blocks, theta, n0=16):
    """在整张单荷环面上证明 Z47 的平均 > theta。"""
    lip = sum(b["lip"] for b in blocks)/4
    pi = arb.pi()
    stack = [(2*pi*(i+arb(1)/2)/n0, 2*pi*(j+arb(1)/2)/n0, pi/n0) for i in range(n0) for j in range(n0)]
    cells = 0
    while stack:
        X, Y, r = stack.pop()
        cells += 1
        val = zone_mean(block_eigs(blocks, X, Y), "Z47")
        if val-lip*r-SLACK > theta:
            continue
        assert val > theta and r > arb("1e-4")
        h = r/2
        stack += [(X+sx*h, Y+sy*h, h) for sx in (-1, 1) for sy in (-1, 1)]
    return cells


LIGHTS = {1: [8, 9, 10], 2: [1, 8, 9], 3: [8, 9, 10], 4: [1, 8, 9], 6: [1, 8, 9]}


def select_rank(d1, d2, d3):
    """Phi(1),Phi(3),Phi(7),Phi(15) = Phi(1)+(0, P1, P2, P3)，P 为增量部分和；唯一严格极小，否则 None。"""
    partial = {1: 0, 3: d1, 7: d1+d2, 15: d1+d2+d3}
    low = min(partial.values())
    winners = [r for r, v in partial.items() if v == low]
    return winners[0] if len(winners) == 1 else None


def zmin_lo_q(lo, q, name):
    return min((lo[(q,)+lab] for lab in ZONES[name]), key=lambda a: a.lower())


def zmax_hi_q(hi, q, name):
    return max((hi[(q,)+lab] for lab in ZONES[name]), key=lambda a: a.upper())


def native_scales(delta):
    """由原壳边表算出全部站点度、站点方差、迹与它们的次序统计量（含软辐边修正）。"""
    geo, soft = shell["geo"], shell["soft"]
    pole, tuned = soft["pole"], soft["star"]["orbit"]
    deg, soft_count = [0]*19, [0]*19
    for i, j, d, z in geo["edges"]:
        for a in (i, j):
            deg[a] += 1
        if pole in (i, j) and (i in tuned or j in tuned):
            soft_count[i] += 1
            soft_count[j] += 1
    assert len(geo["edges"]) == 135
    from collections import Counter
    profile = Counter(zip(deg, soft_count))
    assert profile == {(10, 0): 3, (10, 1): 3, (14, 0): 3, (16, 0): 6, (18, 0): 3, (18, 3): 1}
    weighted = [d-k+k*(1-delta) for d, k in zip(deg, soft_count)]
    variance = [d-k+k*(1-delta)**2 for d, k in zip(deg, soft_count)]
    trace = sum(weighted)
    ordered = sorted(deg)
    assert ordered[9] == 16
    scales = {}
    for (d, k) in sorted(profile):
        scales[f"站点度 ({d},{k}) = {d}-{k}delta"] = d-k*delta
        scales[f"站点方差 ({d},{k}) = {d}-{k}+{k}(1-delta)^2"] = d-k+k*(1-delta)**2
    scales["最小割 T = min 站点方差"] = min(variance, key=lambda a: a.mid())
    scales["最小度"] = min(weighted, key=lambda a: a.mid())
    scales["最大度"] = max(weighted, key=lambda a: a.mid())
    scales["最大站点方差"] = max(variance, key=lambda a: a.mid())
    scales["中位度（无权）16"] = arb(ordered[9])
    scales["中位度（带权）"] = sorted(weighted, key=lambda a: a.mid())[9]
    scales["中位站点方差"] = sorted(variance, key=lambda a: a.mid())[9]
    scales["平均度 tr L0/19"] = trace/19
    scales["平均站点方差"] = sum(variance)/19
    scales["谱中心 c = tr L0/18"] = trace/18
    scales["固定数 9"] = arb(9)
    return scales, trace


def zone_bounds(lo, hi, q, name, upper):
    vals = [(hi if upper else lo)[(q,)+lab] for lab in ZONES[name]]
    return sum(vals)/len(vals)


def verify():
    start = time.time()
    ctx.prec = 128
    basis_list = bases()
    blocks = {q: secular_blocks(q, basis_list) for q in CHARGES}
    lo, hi = certify_bands(blocks)
    for q in CHARGES:
        for blk in blocks[q]:
            assert abs(blk["lip"]-2*arb(2).sqrt()) < arb("1e-12")
    # Delta bound for every positive zero (gamma_1 > 14).
    delta_max = arb(23)/25*arb(2)**-20/(14**2+arb(30))
    assert 4*delta_max < SLACK
    T_lo, T_hi = 9+(1-delta_max)**2, arb(10)

    # ---- A: 谱隙宽度与谱隙比 ----
    def zmin_hi(q, name):
        return min((hi[(q,)+lab] for lab in ZONES[name] if (q,)+lab in hi), key=lambda a: a.upper())
    def zmax_lo(q, name):
        return max((lo[(q,)+lab] for lab in ZONES[name]), key=lambda a: a.lower())
    def zmin_lo(name):
        return min((lo[(q,)+lab] for q in CHARGES for lab in ZONES[name]), key=lambda a: a.lower())
    def zmax_hi(name):
        return max((hi[(q,)+lab] for q in CHARGES for lab in ZONES[name]), key=lambda a: a.upper())
    # Common version: g_r(p) = min_q lambda_{q,r+1}(p) - max_q lambda_{q,r}(p).
    g1 = zmin_lo("Z23")-zmax_hi("Z1")
    rho1 = zmin_lo("Z23")/zmax_hi("Z1")
    pairs = {3: ("Z23", "Z47"), 7: ("Z47", "Z815"), 15: ("Z815", "Z16")}
    gmax = {r: min((zmin_hi(q, b)-zmax_lo(q, a) for q in CHARGES), key=lambda x: x.upper()) for r, (a, b) in pairs.items()}
    rmax = {r: min((zmin_hi(q, b)/zmax_lo(q, a) for q in CHARGES), key=lambda x: x.upper()) for r, (a, b) in pairs.items()}
    assert g1 > arb("4.10589") and rho1 > arb("2.3473")
    assert gmax[3] < arb("3.27791") and gmax[7] < arb("3.93091") and gmax[15] < arb("1.63141")
    assert max(rmax.values(), key=lambda x: x.upper()) < arb("1.4253")
    assert all(g1 > g for g in gmax.values()) and all(rho1 > r for r in rmax.values())
    print("A 共同谱隙（各荷取最小）：秩1 >", g1.lower().str(6), "；秩3,7,15 <", [gmax[r].upper().str(6) for r in (3, 7, 15)])
    print("A 共同谱隙比（各荷取最小）：秩1 >", rho1.lower().str(5), "；秩3,7,15 <", [rmax[r].upper().str(5) for r in (3, 7, 15)])
    # Single-charge version: g_{q,r}(p) = lambda_{q,r+1}(p) - lambda_{q,r}(p).
    single = {}
    for q in CHARGES:
        sg1 = zmin_lo_q(lo, q, "Z23")-zmax_hi_q(hi, q, "Z1")
        sr1 = zmin_lo_q(lo, q, "Z23")/zmax_hi_q(hi, q, "Z1")
        sg = max((zmin_hi(q, b)-zmax_lo(q, a) for a, b in pairs.values()), key=lambda x: x.upper())
        sr = max((zmin_hi(q, b)/zmax_lo(q, a) for a, b in pairs.values()), key=lambda x: x.upper())
        assert sg1 > sg and sr1 > sr
        single[q] = (sg1, sg, sr1, sr)
    SINGLE = {1: ("4.3760", "3.99081", "2.4617", "1.4750"), 2: ("5.4379", "4.13801", "3.6650", "1.4318"),
              3: ("4.2039", "3.94551", "2.3796", "1.4662"), 4: ("5.1926", "4.03601", "3.4362", "1.4253"),
              6: ("4.8730", "3.93091", "3.1373", "1.4280")}
    for q, (a, b, c, d) in SINGLE.items():
        sg1, sg, sr1, sr = single[q]
        assert sg1 > arb(a) and sg < arb(b) and sr1 > arb(c) and sr < arb(d), q
    print("A 单荷谱隙与谱隙比（秩1下界、秩3/7/15上界）：",
          {q: tuple(x.str(5) for x in v) for q, v in single.items()})

    # ---- B: 包含单调泛函的输入 ----
    assert zmin_lo("Z1") > arb("0.0411") and zmin_lo("Z23") > arb("7.1530")
    print("B 全部本征值 >", zmin_lo("Z1").lower().str(4), "，第二层起 >", zmin_lo("Z23").lower().str(5))

    # ---- C: 费米能级狄拉克海 E_mu(r)=sum_{j<=r}(lambda_j-mu) ----
    bounds = {}
    for name in ("Z23", "Z47", "Z815"):
        bounds[name] = (min((zone_bounds(lo, hi, q, name, False) for q in CHARGES), key=lambda a: a.lower()),
                        max((zone_bounds(lo, hi, q, name, True) for q in CHARGES), key=lambda a: a.upper()))
    tight = arb("10.6")
    cells = sum(adaptive_inf_a47(blocks[q], tight) for q in CHARGES)
    inf47 = max(bounds["Z47"][0], tight, key=lambda a: a.lower())
    W1 = bounds["Z23"][0]
    W3 = (bounds["Z23"][1], inf47)
    W7 = (bounds["Z47"][1], bounds["Z815"][0])
    W15 = bounds["Z815"][1]
    assert W1 > arb("7.37474") and W3[0] < arb("8.22501") and W3[1] > arb("10.5999")
    assert cells == 35752
    assert T_lo-W3[0] > arb("1.774") and W3[1]-T_hi > arb("0.5999")
    assert W7[0] < arb("11.22305") and W7[1] > arb("15.81795") and W15 < arb("16.59516")
    print("C 一致窗：秩1 mu<", W1.lower().str(6), " 秩3", W3[0].upper().str(6), W3[1].lower().str(5),
          " 秩7", W7[0].upper().str(6), W7[1].lower().str(6), " 秩15 mu>", W15.upper().str(6), " 自适应格数", cells)
    # Witnesses on the single-charge torus; each extends to an open set.
    pi = arb.pi()
    wit = {}
    for key, q, a, b, name in (("inf23", 3, 128, 128, "Z23"), ("sup23", 2, 247, 9, "Z23"),
                               ("inf47", 6, 96, 158, "Z47"), ("sup47", 1, 0, 0, "Z47"),
                               ("inf815", 3, 0, 0, "Z815"), ("sup815", 4, 0, 0, "Z815"),
                               ("inf815b", 2, 128, 128, "Z815")):
        wit[key] = zone_mean(block_eigs(blocks[q], 2*pi*a/256, 2*pi*b/256), name)
    assert wit["inf23"] < arb("7.37523") and wit["sup23"] > arb("8.22064")
    assert wit["inf47"] < arb("10.63563") and wit["sup47"] > arb("11.22211")
    assert wit["inf815"] < arb("15.96474") and wit["sup815"] > arb("16.59361") and wit["inf815b"] > arb("16.47826")
    print("C 见证点：", {k: v.str(8) for k, v in wit.items()})

    # Native scales.
    # Native scales computed from the edge list; delta ranges over [0, delta_max].
    delta_ball = arb(delta_max/2, delta_max/2)
    scales, trace = native_scales(delta_ball)
    for dv in (arb(0), delta_max):
        sc, tr = native_scales(dv)
        assert abs(tr-(270-6*dv)) < arb("1e-30")
        assert abs(sc["谱中心 c = tr L0/18"]-(15-dv/3)) < arb("1e-30")
        assert abs(sc["最小割 T = min 站点方差"]-(9+(1-dv)**2)) < arb("1e-30")
    assert abs(scales["平均度 tr L0/19"]-arb("14.2105")) < arb("1e-4")
    def rank_of(x):
        if x < W1: return 1
        if x > W3[0] and x < W3[1]: return 3
        if x > W7[0] and x < W7[1]: return 7
        if x > W15: return 15
        if x > wit["inf815"] and x < wit["inf815b"]: return "7/15"
        return None
    chosen = {k: rank_of(v) for k, v in scales.items()}
    expected = {}
    for k in scales:
        if k.startswith(("站点度 (10", "站点方差 (10", "最小割", "最小度", "固定")):
            expected[k] = 3
        elif k.startswith(("站点度 (14", "站点方差 (14", "平均", "谱中心")):
            expected[k] = 7
        elif k.startswith(("站点度 (16", "站点方差 (16", "中位")):
            expected[k] = "7/15"
        else:
            expected[k] = 15
    assert chosen == expected, chosen
    # Scale 16 lies strictly between two witness values of a815: both ranks win on open sets.
    assert wit["inf815"] < 16 < wit["inf815b"]
    print("C 原生能级的选择：")
    for k, v in scales.items():
        print("   ", k, v.str(10), "->", chosen[k])

    # psi-generalization at mu=T: power means M_s with s >= -20, Boltzmann beta <= 2.
    s, beta = -20, 2
    for q in CHARGES:
        low47 = [lo[(q,)+lab] for lab in ZONES["Z47"]]
        assert (sum(x**s for x in low47)/4)**(arb(1)/s) > T_hi
        assert sum((-beta*(x-T_hi)).exp() for x in low47)/4 < 1
    wq = block_eigs(blocks[2], pi, pi)
    w47 = [wq[k][i-1] for k, i in ZONES["Z47"]]
    lam4 = min(w47, key=lambda a: a.mid())
    assert arb("9.60926") < lam4 < arb("9.60928") and lam4 < T_lo
    assert sum(1 for e in wq for x in e if x < T_lo) == 4
    assert sum((-4*(x-T_lo)).exp() for x in w47)/4 > 1
    assert (sum(x**-40 for x in w47)/4)**(arb(-1)/40) < T_lo
    assert wit["sup23"] < T_lo and zone_mean(wq, "Z47") > T_hi
    print("C 推广：mu=T 下幂平均 s>=-20 与玻尔兹曼 beta<=2 处处选3；(pi,pi) 荷二处 lambda4=", lam4.str(8),
          "，beta=4 与 s=-40 失败；逐模 T 截断在此给 4 个轻模")

    # Actual zeros: full matrices with soft spokes.
    heights = acb.zeta_zeros(1, 30)
    lap = shell["laplacian"]
    per_mode_fail = []
    for n, z in enumerate(heights, start=1):
        delta = arb(23)/25*arb(2)**-20/(z.imag**2+shell["soft"]["point"]["study"]["offset"])
        T = 9+(1-delta)**2
        for q in CHARGES:
            ev = sorted((v.real for v in lap(q, z.imag).eig()), key=lambda a: a.mid())
            m23, m47, m815 = (sum(ev[a:b])/(b-a) for a, b in ((1, 3), (3, 7), (7, 15)))
            assert m23 < T < m47 and m47 < 15-delta/3 < m815
            assert W3[0] > m23 and m47 > W3[1]-SLACK
            below = sum(1 for v in ev if v < T)
            assert all(v < T or v > T for v in ev)
            if below != 3:
                per_mode_fail.append((n, q, below))
    assert per_mode_fail == [(3, 4, 4), (10, 2, 4), (12, 2, 4), (13, 6, 4), (20, 6, 4)]
    print("前 30 枚零点：Fermi 能级 T 恒选 3、c 恒选 7；逐模 T 截断失效于", per_mode_fail)
    # Two-step reading h=|lambda(lambda-c)| folds near lambda=0 and lambda=c: at U=V=1.
    lights = {}
    for q in CHARGES:
        ev = sorted((x for e in block_eigs(blocks[q], arb(0), arb(0)) for x in e), key=lambda a: a.mid())
        assert all(abs(x*(x-15)) < T_lo or abs(x*(x-15)) > T_hi for x in ev)
        lights[q] = [j for j, x in enumerate(ev, start=1) if abs(x*(x-15)) < T_lo]
    assert lights == LIGHTS, lights
    print("两步读数 h<T 在 U=V=1 处所取的单步层号", lights)

    # Complete discrimination for Phi(r)=sum_{j<=r} phi(lambda_j): increments over zones.
    # Counterexample: phi=-1 (lambda<9), +1 (9<=lambda<13), -1 (lambda>=13).
    assert arb(ZONE_WINDOWS["Z23"][1]) < 9 and arb(ZONE_WINDOWS["Z47"][0]) >= 9
    assert arb(ZONE_WINDOWS["Z47"][1]) < 13 and arb(ZONE_WINDOWS["Z815"][0]) >= 13
    step = (-2, 4, -8)
    assert select_rank(*step) == 15
    # Block averages of phi change sign from Z23 (-1) to Z47 (+1), yet rank 15 wins everywhere.
    ev = sorted((x for e in block_eigs(blocks[1], arb(0), arb(0)) for x in e), key=lambda a: a.mid())
    phi = lambda x: -1 if x < 9 else (1 if x < 13 else -1)
    totals = [sum(phi(x) for x in ev[:r]) for r in (1, 3, 7, 15)]
    assert totals == [-1, -3, 1, -7]
    # phi(lambda)=lambda-c: Delta_1<0, Delta_2<0, Delta_3>0 everywhere selects 7.
    assert W7[0] < 15-delta_max and 15 < W7[1]
    assert select_rank(-1, -1, 1) == 7 and select_rank(-1, 1, 1) == 3 and select_rank(1, 1, 1) == 1
    assert select_rank(-1, 1, -1) is None and select_rank(-1, 1, -2) == 15
    print("可加泛函判别：阶梯 phi 的 Phi(1,3,7,15) =", totals, "，选 15")
    print("用时 %.1f 秒" % (time.time()-start))


if __name__ == "__main__":
    verify()
