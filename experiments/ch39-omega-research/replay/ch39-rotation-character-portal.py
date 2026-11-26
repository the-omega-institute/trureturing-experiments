"""由原壳三阶轮换生成的群代数门户能否在四枚规范候选上抵消轻 Higgs 二次项的两圈有限阈值。

门户取为轮换群代数中的厄米元 P = a I + b S + conj(b) S^2，作用在荷三重标量双重态上。
它在特征标 chi 的本征线上的对角元恰为 p_chi，三元实数 (p_0,p_1,p_2) 与 (a,b) 一一对应。
一圈门户贡献为 sum_chi p_chi W_chi，其中 W_chi 是该特征标的活跃重双重态的
m^2 (1 + 2 log(Q/m)) / (16 pi^2) 之和；抵消要求它等于纯电弱两圈有限阈值 T。
"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

ROOT = Path(__file__).parent
common = run_path(str(ROOT / "ch39-chiral-common-portal-scale.py"))
fixed, chiral, two = common["fixed"], common["chiral"], common["two"]
geometry = run_path(str(ROOT / "ch39-fixed-planck-magnetic-window.py"))["geometry"]
matter = run_path(str(ROOT / "ch39-charge-matter-determinant-order.py"))
soft = matter["soft"]

ZEROS = (11661, 11934, 12207, 28325)
SCALE = arb("3.725e14")
UNIT = arb("2.435e18")/724**2
TOLERANCE = arb("125.25")**2/2
ALPHA = {11661: ".11838", 11934: ".11846", 12207: ".11837", 28325: ".11838"}
# 说明文档表中的四舍五入值：T/1e26 与 W/1e27。
TABLE = {11661: ("1.04651", ("2.613824", "2.760040", "1.833693")),
         11934: ("1.04518", ("2.608613", "2.736231", "1.818032")),
         12207: ("1.04669", ("2.614174", "2.764174", "1.835853")),
         28325: ("1.04661", ("2.614720", "2.759766", "1.835275"))}


def rotation():
    """坐标轮换 (a,b,c)->(b,c,a) 的站点置换及其三个特征标子空间的标准正交基。"""
    sites = geometry["sites"]
    perm = [sites.index((b, c, a)) for a, b, c in sites]
    assert all(perm[perm[perm[i]]] == i for i in range(19))
    assert [i for i in range(19) if perm[i] == i] == [sites.index((17, 17, 17))]
    edges = {(i, j): (d, u) for i, j, d, u in geometry["edges"]}
    for (i, j), (d, u) in edges.items():
        ri, rj = perm[i], perm[j]
        # 有钟边保持方向，因而同一钟相位落在同向边上。
        if d in (724, 804):
            assert ri < rj
        assert edges[tuple(sorted((ri, rj)))][0] == d
        # 整数核对：内积与以极点为首行的带号行列式不变，故原生相位不变。
        pole = (17, 17, 17)
        assert geometry["dot"](sites[ri], sites[rj]) == geometry["dot"](sites[i], sites[j])
        assert geometry["det"](pole, sites[ri], sites[rj]) == geometry["det"](pole, sites[i], sites[j])
    S = acb_mat([[int(i == perm[j]) for j in range(19)] for i in range(19)])
    omega = acb(-arb(1)/2, arb(3).sqrt()/2)
    orbits, seen = [], set()
    for i in range(19):
        if i not in seen:
            orbit = [i] if perm[i] == i else [i, perm[i], perm[perm[i]]]
            seen.update(orbit)
            orbits.append(orbit)
    assert sorted(len(o) for o in orbits) == [1, 3, 3, 3, 3, 3, 3]
    bases = []
    for chi in range(3):
        cols = [o for o in orbits if len(o) == 3 or chi == 0]
        B = acb_mat(19, len(cols))
        for c, o in enumerate(cols):
            for k, i in enumerate(o):
                B[i, c] = omega**(-chi*k)/arb(len(o)).sqrt()
        assert all(abs(z) < arb("1e-50") for z in (S*B-B*omega**chi).entries())
        bases.append(B)
    assert [b.ncols() for b in bases] == [7, 6, 6]
    return S, bases


def charge_matrix(q, height):
    # 与 spectral_supplier 同一十二条有钟边的单位电导矩阵。
    clocks = geometry["phases"](height)
    mat = acb_mat(19, 19)
    for i, j, d, u in geometry["edges"]:
        mat[i, i] += 1
        mat[j, j] += 1
        value = -u**q
        if d in (724, 804):
            value = value*clocks[d]**q
        mat[i, j] = value
        mat[j, i] = value.conjugate()
    return mat


def spoke_removal(S, height):
    """三条辐边（极点到一个轮换轨道）的刚度减量，每边权重 amount，传输为一。"""
    pole, orbit = soft["pole"], soft["star"]["orbit"]
    assert geometry["sites"][pole] == (17, 17, 17)
    star = acb_mat(19, 19)
    for v in orbit:
        e = [int(i == v)-int(i == pole) for i in range(19)]
        for i in range(19):
            for j in range(19):
                star[i, j] += e[i]*e[j]
    assert star.trace().real == 6
    perm = [geometry["sites"].index((b, c, a)) for a, b, c in geometry["sites"]]
    assert perm[pole] == pole and {perm[v] for v in orbit} == set(orbit) and len(set(orbit)) == 3
    amount = arb(23)/25*arb(2)**-20/(height**2+soft["point"]["study"]["offset"])
    return amount, star


def candidate(model, initial, read, S, bases, n):
    height = acb.zeta_zero(n).imag
    rows = read(height)
    levels = {q: chiral["ordered"]([h for h, z in rows if z == q]) for q in (0, *model["charges"])}
    assert all(levels[q][2] < 9 < levels[q][3] for q in model["charges"])
    _, _, alpha = fixed["crossing_enclosure"](model, initial, levels, False, 724, "29.040", ".1184")
    assert arb(".117") < alpha < arb(".119")
    table = [((UNIT/model["mz"]).log()+l, w, B, f, s) for l, w, B, f, s
             in two["thresholds"](model, levels, False)]
    target, moments = common["portal_at"](model, initial, (levels, 1/alpha, table), SCALE, steps=65536)
    # 交点模块会改动工作精度，显式算子在 192 位下重算。
    prec, ctx.prec = ctx.prec, 192
    unit_matrix = charge_matrix(3, height)
    # 与 shell_spectra 的分解一致：极点边单位传输的星图加非辐边荷三矩阵。
    base, neutral, _ = soft["vector_pencil"](height)
    assert all(z.contains(0) for z in (unit_matrix-base-matter["nonspokes"](3, height)).entries())
    amount, star = spoke_removal(S, height)
    mat = unit_matrix-amount*star
    # 含辐边减量的荷三算子与 S 对易：每个矩阵元的区间都包含零。
    assert all(z.contains(0) for z in (S*mat-mat*S).entries())
    # 中性算子的迹为 270，减量迹为 6a，且常向量为零模，故十八个正本征值的均值为 15-a/3。
    trace = (base+neutral).trace().real
    assert trace == 270
    center = (trace-6*amount)/18
    labelled = []
    for chi, B in enumerate(bases):
        for eig in (B.conjugate().transpose()*mat*B).eig(multiple=True, algorithm="rump"):
            assert eig.imag.contains(0)
            v = eig.real
            labelled.append((v*abs(v-center), chi, v))
    assert len(labelled) == 19
    gap = min(abs(a[2]-b[2]) for a in labelled for b in labelled if a[1] != b[1])
    labelled.sort(key=lambda t: float(t[0].mid()))
    assert all(a[0] < b[0] for a, b in zip(labelled, labelled[1:]))
    # 显式算子的包围落在 spectral_supplier 的 1e-12 宽包围之内，次序一致。
    assert all(g.overlaps(h) for (g, _, _), h in zip(labelled, levels[3]))
    W = [arb(0)]*3
    x = [arb(0)]*4
    pattern, heights = [], []
    for g, chi, _ in labelled[3:]:
        mass = UNIT*g
        assert mass < SCALE or mass > SCALE
        if mass < SCALE:
            w = mass**2*(1+2*(SCALE/mass).log())/(16*arb.pi()**2)
            W[chi] += w
            for k in range(4):
                x[k] += w*(mass/(80*UNIT))**(2*k)
            pattern.append(chi)
            heights.append(g)
    assert len(pattern) == 14
    ctx.prec = prec
    assert all((a-b).contains(0) or a.overlaps(b) for a, b in zip(x, moments[:4]))
    return alpha, target, W, x, pattern, heights, gap


def cone_interval(Wa, Wb):
    ratios = [b/a for a, b in zip(Wa, Wb)]
    lo, hi = ratios[0], ratios[0]
    for r in ratios[1:]:
        lo, hi = lo.min(r), hi.max(r)
    return lo, hi


def null_residual(rows, targets):
    """4x3 方程的左零向量 y；min_p max_n |T_n-(Wp)_n| 恰为 |y.T|/||y||_1。"""
    y = []
    for i in range(4):
        minor = arb_mat([r for k, r in enumerate(rows) if k != i])
        y.append((-1)**i*minor.det())
    for j in range(3):
        assert sum(y[i]*rows[i][j] for i in range(4)).contains(0)
    return abs(sum(a*t for a, t in zip(y, targets)))/sum(abs(a) for a in y), y


def psd_lower_bound(rows, targets):
    """min_{p>=0} max_n |T_n-(Wp)_n| 的对偶下界。"""
    import numpy as np
    from scipy.optimize import linprog
    k = len(rows)
    W = np.array([[float(w.mid())/1e27 for w in r] for r in rows])
    T = np.array([float(t.mid())/1e26 for t in targets])
    c = np.r_[-T, T]
    A = np.vstack([np.r_[np.ones(k), np.ones(k)], np.hstack([W.T, -W.T])])
    res = linprog(c, A_ub=A, b_ub=np.r_[1, np.zeros(3)], bounds=[(0, None)]*(2*k))
    assert res.status == 0
    y = [arb(float(a))-arb(float(b)) for a, b in zip(res.x[:k], res.x[k:])]
    # 对任意 p>=0：y.T <= ||y||_1 e + sum_chi max(0,(W^T y)_chi) p_chi。
    # 若 e <= min T，则 (Wp)_n <= 2 T_n，于是 p_chi <= 2 T_n / W_{n chi}。
    slack = arb(0)
    for chi in range(3):
        g = sum(y[n]*rows[n][chi] for n in range(k))
        pmax = min((2*targets[n]/rows[n][chi]).upper() for n in range(k))
        slack += arb(g.upper()).max(0)*pmax
    norm = sum(abs(a) for a in y)
    value = (sum(a*t for a, t in zip(y, targets))-slack)/norm
    return value.min(min(targets))


def verify():
    model, _, _, initial, _ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    S, bases = rotation()
    data, patterns, clusters = {}, {}, {}
    for n in ZEROS:
        alpha, T, W, x, pattern, heights, gap = candidate(model, initial, read, S, bases, n)
        clusters[n] = heights
        print("不同特征标本征值的最小间距", n, gap.lower().str(6), flush=True)
        if n == 11934:
            assert arb("3.45e-4") < gap < arb("3.47e-4")
        else:
            assert arb(".015") < gap < arb(".018")
        assert abs(alpha-arb(ALPHA[n])) < arb("1e-5")
        assert [pattern.count(c) for c in range(3)] == [5, 5, 4]
        assert all(w > 0 for w in W) and T > 0 and arb("8.5e19") < T.rad() < arb("9.7e19")
        ref_T, ref_W = TABLE[n]
        assert abs(T-arb(ref_T)*10**26) < arb("5e20")
        assert all(abs(w-arb(r)*10**27) < arb("1e21") for w, r in zip(W, ref_W))
        data[n] = (T, W, x)
        patterns[n] = pattern
        print("零点、强耦合、两圈阈值 T/GeV²、三特征标权 W/GeV²、活跃层特征标",
              n, alpha.str(8), T.str(10), [w.str(9) for w in W], pattern, flush=True)

    # 不同特征标块之间无能级排斥：11934 在 h 约 37 与 44 的两簇处交换了特征标次序。
    assert patterns[11661] == patterns[12207] == patterns[28325]
    diff = [i for i in range(14) if patterns[11934][i] != patterns[11661][i]]
    assert diff == [4, 5, 7, 8, 9]
    for n in ZEROS:
        assert all(arb(36) < clusters[n][i] < arb(38) for i in (4, 5))
        assert all(arb(44) < clusters[n][i] < arb(45) for i in (7, 8, 9))

    # 一、两枚、非负特征标门户：可行当且仅当 T_b/T_a 落在三特征标权比的闭区间内。
    feasible = {}
    for a, b in combinations(ZEROS, 2):
        lo, hi = cone_interval(data[a][1], data[b][1])
        r = data[b][0]/data[a][0]
        if lo < r < hi:
            feasible[(a, b)] = True
        else:
            assert r > hi or r < lo
            feasible[(a, b)] = False
            # 非负 p 时 (Wb.p)/(Wa.p) 在 [lo,hi] 内，最佳残差下界取比值间隙。
            gap = (r-hi) if r > hi else (lo-r)
            residual = gap*data[a][0]/(1+(hi if r > hi else lo))
            assert residual > {(11661, 11934): arb("3.753e22"), (11934, 12207): arb("3.544e22"),
                               (11934, 28325): arb("5.057e22")}[(a, b)]
            print("非负特征标门户排除", a, b, "比值、区间、两枚最大残差下界/GeV²",
                  r.str(10), lo.str(10), hi.str(10), residual.lower().str(6))
    assert feasible == {(11661, 11934): False, (11661, 12207): True, (11661, 28325): True,
                        (11934, 12207): False, (11934, 28325): False, (12207, 28325): True}
    print("两枚非负特征标门户的可行性", feasible)

    # 二、三枚：唯一解，第二个非平凡特征标的系数恒为负；目标放宽电弱容差后仍成立。
    def widened(n):
        return data[n][0]+arb(0, TOLERANCE.upper())
    for triple in combinations(ZEROS, 3):
        M = arb_mat([data[n][1] for n in triple])
        p = M.solve(arb_mat([[widened(n)] for n in triple])).entries()
        assert p[0] > 0 and p[1] > 0 and p[2] < 0
        assert arb(".034") < p[0] < arb(".046") and arb(".007") < p[1] < arb(".020")
        assert arb("-.036") < p[2] < arb("-.015")
        print("三枚的唯一特征标门户 (p_0,p_1,p_2)", triple, [v.str(5) for v in p])

    # 三、半正定特征标门户 (p>=0) 的极小极大残差：浮点线性规划给对偶向量，Arb 复核。
    for chosen in [*combinations(ZEROS, 3), ZEROS]:
        bound = psd_lower_bound([data[n][1] for n in chosen], [data[n][0] for n in chosen])
        assert bound > (arb("9.25e21") if chosen == (11661, 12207, 28325) else arb("3.75e22"))
        if 11934 in chosen and 28325 in chosen:
            assert bound > arb("5.057e22")
        print("半正定特征标门户的极小极大残差下界/GeV²", chosen, bound.lower().str(6))

    lows = []
    # 三之二、实门户 b 为实数，等价 p_1=p_2：二维族在任意三枚上被排除。
    for triple in combinations(ZEROS, 3):
        rows = [[data[n][1][0], data[n][1][1]+data[n][1][2]] for n in triple]
        yy = [(-1)**i*arb_mat([r for k, r in enumerate(rows) if k != i]).det() for i in range(3)]
        low = abs(sum(a*data[n][0] for a, n in zip(yy, triple)))/sum(abs(a) for a in yy)
        assert low > arb("7.614e21") and arb(low.lower()) < arb("2.1e22")
        lows.append(low)
        print("实轮换门户 (p_1=p_2) 在三枚上的极小极大残差下界/GeV²", triple, low.lower().str(6))

    assert max(v.lower() for v in lows) > arb("2.092e22")

    # 四、四枚：有符号特征标门户被排除，残差为精确极小极大值 |y.T|/||y||_1。
    targets = [data[n][0] for n in ZEROS]
    char_res, y = null_residual([data[n][1] for n in ZEROS], targets)
    assert char_res > arb("3.191e21")
    assert char_res/TOLERANCE > arb("4e17")
    print("四枚共用有符号特征标门户的极小极大残差下界/GeV²", char_res.lower().str(6),
          "左零向量", [v.str(5) for v in y])

    # 四之二、显式辐边减量使矩的包围收紧，四枚上的有符号二次谱门户也可判定为排除；
    # 四阶矩阵非奇异，三次谱门户在四枚上唯一可解。
    qres, _ = null_residual([data[n][2][:3] for n in ZEROS], targets)
    assert qres > arb("1.787e21")
    det = arb_mat([data[n][2][:4] for n in ZEROS]).det()
    assert not det.contains(0)
    print("四枚共用二次谱门户的极小极大残差下界/GeV²、三次矩阵行列式", qres.lower().str(6), det.str(5))
    # 三次门户的系数包围（右端取 T_n±ε）；它在四枚全部活跃层上取正值，但首项为负，实轴上不处处为正。
    # 用克莱姆法则：通用求解器在此近共线矩阵上给出过宽的包围。
    rhs = [data[n][0]+arb(0, TOLERANCE.upper()) for n in ZEROS]
    # 按右端所在列展开，余子式只含矩；门户值 p(z) 对 T 线性，逐层以紧的系数计算，避免依赖放大。
    M = [data[n][2][:4] for n in ZEROS]
    cof = [[(-1)**(i+k)*arb_mat([[M[r][j] for j in range(4) if j != k]
                                 for r in range(4) if r != i]).det()/det for k in range(4)]
           for i in range(4)]
    cubic = [sum(rhs[i]*cof[i][k] for i in range(4)) for k in range(4)]
    print("四枚的唯一三次门户系数上下界", [(c.lower().str(5), c.upper().str(5)) for c in cubic])
    bounds = (("0.0399", "0.0410"), ("-0.1536", "-0.1464"), ("0.2607", "0.2800"), ("-0.1677", "-0.1528"))
    assert all(arb(lo) < c < arb(hi) for c, (lo, hi) in zip(cubic, bounds))
    least = None
    for n in ZEROS:
        for g in clusters[n]:
            z = (g/80)**2
            value = sum(rhs[i]*sum(cof[i][k]*z**k for k in range(4)) for i in range(4))
            least = value if least is None else least.min(value)
    assert least > arb(".0037")
    assert cubic[3] < 0
    print("三次门户在四枚 56 个活跃层上的最小值下界", least.lower().str(5))

    # 五、有符号谱多项式：沿 11661,11934,12207 依次加入，次数 N-2 被排除，次数 N-1 唯一可解。
    for N in (2, 3):
        chosen = ZEROS[:N]
        rows = [data[n][2][:N-1] for n in chosen]
        yy = [(-1)**i*arb_mat([r for k, r in enumerate(rows) if k != i]).det() for i in range(N)]
        low = abs(sum(a*data[n][0] for a, n in zip(yy, chosen)))/sum(abs(a) for a in yy)
        assert low > (arb("2.583e23") if N == 2 else arb("1.418e22"))
        full = arb_mat([data[n][2][:N] for n in chosen])
        coef = full.solve(arb_mat([[widened(n)] for n in chosen])).entries()
        print("零点数、次数 N-2 的极小极大残差下界/GeV²、次数 N-1 的唯一系数",
              N, low.lower().str(6), [c.str(5) for c in coef])
    # 三枚的二次门户：系数与已知结果一致，且在整个实轴上为正（半正定谱门户）。
    c0, c1, c2 = coef
    assert arb(".022") < c0 < arb(".024") and arb("-.066") < c1 < arb("-.062") and arb(".067") < c2 < arb(".073")
    assert c0-c1**2/(4*c2) > arb(".0077")
    print("二次门户实轴极小值下界", (c0-c1**2/(4*c2)).lower().str(5))

if __name__ == "__main__":
    verify()
