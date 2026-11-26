"""两步质量读数下容许秩含 2、3、4 时五荷一致的三代费米能级窗为空，而全环面共同两步分界只剩秩 15。

说明：两步质量读数下容许秩含二三四时五荷一致的三代费米能级窗为空而全环面共同两步分界只剩秩15.md

模型：十九点贝伦德球壳，五种绝对荷 q in {1,2,3,4,6}，两钟相位 (U,V) 取遍整张环面，
辐边亏额 0<=delta<delta_max=(23/25)2^-20/(14^2+30)。第七章定义 7.1 的两步质量算子
D_q=L_q(L_q-cI)，c=tr L_0/18=15-delta/3，读数 h_{q,j}=|lambda_j(lambda_j-c)|，每个荷内按 h 升序
编号 h_(1)<=...<=h_(19)（第四章 4.2 节“按 h 升序切出 3 个”的口径）。
相位约定：build_family 的 (U,V) 是荷 q 自己的相位 (u^q,v^q)。凡“五荷在同一点”的陈述，
点指联合相位 (u,v)=(e^{iX},e^{iY})，荷 q 代入 (qX,qY)；只对单个荷取全环面上下确界的陈述
不受影响，因为 (u,v)->(u^q,v^q) 是环面到自身的满射。

证书。
A. 标签带界（复用 ch39-native-rank-selector-classes.py 的 Haynsworth 全环面证书）给出每个
   单步标签在整张环面上的区间；经 h 的精确值域（端点、临界点 c/2、零点 0 与 c）变为 h 区间。
   若 x_i in [l_i,u_i]，则第 j 个次序统计量落在 [l_(j),u_(j)]。由此证明：秩 15 是整张环面上
   五荷逐点共同的两步分界；荷一、荷三各自有一致的三模窗。
B. 见证点：在 X,Y in 2pi Z/64 或 2pi Z/4096 处用 Arb 求本征值。严格不等式在邻域内保持。
   (i) 秩 r in {1..18}\\{15} 在某个联合点上不是五荷共同分界；(ii) 荷二、四、六的 h_(3) 在一点
   大于 h_(4) 在另一点的值，由环面连通与介值定理，两者值域相交。
C. 自适应网格：格心 Arb 本征值，块内扰动范数 <= (|b_1|+|b_2|) r（Weyl），软辐边下移 <=4 delta。
   在单步容许秩 {1,3,7,15} 下证明 sup a23 < 13、inf a47 > 17（a 为 h 次序统计量的层平均）。
D. 前 100 枚实际正零点（含各自的 delta）逐点 Arb 检验。
"""
from pathlib import Path
from runpy import run_path
import time

from flint import acb, arb, ctx

HERE = Path(__file__).with_name
prev = run_path(str(HERE("ch39-native-rank-selector-classes.py")))
shell = prev["shell"]
CHARGES = (1, 2, 3, 4, 6)
DMAX = arb(23)/25*arb(2)**-20/(14**2+arb(30))
LAMSL = arb("1e-7")          # 4 delta_max < 1e-7
GERSH = arb(36)              # lambda <= 2 * 最大度 18


def emin(xs):
    m = xs[0]
    for x in xs[1:]:
        if x < m:
            m = x
    return m


def emax(xs):
    m = xs[0]
    for x in xs[1:]:
        if x > m:
            m = x
    return m


def hrange(a, b, C):
    """lambda in [a,b]（精确端点）、c in C 时 |lambda(lambda-c)| 的包围 [lo,hi]，端点精确。"""
    cand = [abs(a*(a-C)), abs(b*(b-C))]
    if not (b < C/2) and not (a > C/2):
        cand.append(C*C/4)
    if not (b < C) and not (a > C):
        cand.append(arb(0))
    if not (b < 0) and not (a > 0):
        cand.append(arb(0))
    return emin([x.lower() for x in cand]), emax([x.upper() for x in cand])


def order_stats(intervals):
    """x_i in [l_i,u_i] => x_(j) in [l_(j),u_(j)]；返回两列已排序的精确端点。"""
    return sorted(a for a, b in intervals), sorted(b for a, b in intervals)


def block_eigs(q, X, Y):
    U, V = acb(0, X).exp(), acb(0, Y).exp()
    out = []
    for blk in BLOCKS[q]:
        P = U*blk["terms"][0]+V*blk["terms"][1]
        M = blk["block"]+P+P.conjugate().transpose()
        out.append([x.real for x in M.eig(multiple=True)])
    return out


def cell_stats(q, ev, r, C):
    """格心本征值 ev、格半径 r：格内任一点、任一 delta 下 h 次序统计量的包围。"""
    iv = []
    for blk, e in zip(BLOCKS[q], ev):
        eps = blk["lip"]*r
        los = sorted((x-eps-LAMSL).lower() for x in e)
        his = sorted((x+eps).upper() for x in e)
        iv += [hrange(a, b, C) for a, b in zip(los, his)]
    assert len(iv) == 19
    return order_stats(iv)


def point_stats(q, X, Y, C):
    return cell_stats(q, block_eigs(q, X, Y), arb(0), C)


def adaptive(q, A, B, C, n0=16):
    """在整张单荷环面上证明 (h_(2)+h_(3))/2 < A 且 (h_(4)+...+h_(7))/4 > B。"""
    pi = arb.pi()
    stack = [(2*pi*(i+arb(1)/2)/n0, 2*pi*(j+arb(1)/2)/n0, pi/n0, True, True)
             for i in range(n0) for j in range(n0)]
    cells = 0
    while stack:
        X, Y, r, pa, pb = stack.pop()
        cells += 1
        ev = block_eigs(q, X, Y)
        lo, hi = cell_stats(q, ev, r, C)
        if pa and (hi[1]+hi[2])/2 < A:
            pa = False
        if pb and (lo[3]+lo[4]+lo[5]+lo[6])/4 > B:
            pb = False
        if pa or pb:
            lo0, hi0 = cell_stats(q, ev, arb(0), C)
            assert not pa or (lo0[1]+lo0[2])/2 < A, ("a23 bound fails", q)
            assert not pb or (hi0[3]+hi0[4]+hi0[5]+hi0[6])/4 > B, ("a47 bound fails", q)
            assert r > arb("1e-4")
            h = r/2
            stack += [(X+sx*h, Y+sy*h, h, pa, pb) for sx in (-1, 1) for sy in (-1, 1)]
    return cells


def grid(m, n):
    return 2*arb.pi()*m/n


def verify():
    global BLOCKS
    t0 = time.time()
    ctx.prec = 64
    C = (arb(15)-DMAX/3).union(arb(15))     # c = 15 - delta/3
    basis = prev["bases"]()
    BLOCKS = {q: prev["secular_blocks"](q, basis) for q in CHARGES}
    lo_b, hi_b = prev["certify_bands"](BLOCKS)

    # ---- A. 标签带界 -> 两步次序统计量的全环面包围 ----
    band = {}
    for q in CHARGES:
        iv = []
        for k, size in ((0, 7), (1, 6), (2, 6)):
            for i in range(1, size+1):
                a = lo_b[q, k, i].lower()
                b = hi_b[q, k, i].upper() if (q, k, i) in hi_b else GERSH
                iv.append(hrange(a, b, C))
        band[q] = order_stats(iv)
    u15 = emax([band[q][1][14] for q in CHARGES])
    l16 = emin([band[q][0][15] for q in CHARGES])
    assert arb("60.3921") < u15 < arb("60.3922") and arb("61.9199") < l16 < arb("61.9200")
    assert u15 < l16 and u15 < arb("60.39217") and l16 > arb("61.91994")
    # 两步秩 15 取到的恰是单步第 1..15 层的十五个标签：其 h 上界都低于第 16..19 层标签的 h 下界。
    for q in CHARGES:
        iv = {}
        for k, size in ((0, 7), (1, 6), (2, 6)):
            for i in range(1, size+1):
                b = hi_b[q, k, i].upper() if (q, k, i) in hi_b else GERSH
                iv[k, i] = hrange(lo_b[q, k, i].lower(), b, C)
        top = prev["ZONES"]["Z16"]
        assert len(top) == 4
        assert all(v[1] < iv[t][0] for x, v in iv.items() if x not in top for t in top)
    print("A 秩 15：整张环面、五荷逐点 max_q h_(15) <", u15.str(7), "<", l16.str(7), "< min_q h_(16)")
    win13 = {}
    for q in (1, 3):
        u3, l4 = band[q][1][2], band[q][0][3]
        assert u3 < l4
        win13[q] = (u3, l4)
    assert win13[1][0] < arb("0.167877") and win13[1][1] > arb("6.54293")
    assert win13[3][0] < arb("0.909794") and win13[3][1] > arb("7.04633")
    print("A 荷一：sup h_(3) <", win13[1][0].str(6), "，inf h_(4) >", win13[1][1].str(6),
          "；荷三：sup h_(3) <", win13[3][0].str(6), "，inf h_(4) >", win13[3][1].str(6))
    # 第 1..3 层在荷一、荷三上恰为单步标签 (0,4),(1,3),(2,3)：其 h 上界都低于其余标签的 h 下界。
    for q in (1, 3):
        iv = {}
        for k, size in ((0, 7), (1, 6), (2, 6)):
            for i in range(1, size+1):
                b = hi_b[q, k, i].upper() if (q, k, i) in hi_b else GERSH
                iv[k, i] = hrange(lo_b[q, k, i].lower(), b, C)
        light = [(0, 4), (1, 3), (2, 3)]
        assert emax([iv[x][1] for x in light]) < emin([v[0] for x, v in iv.items() if x not in light])

    # ---- B. 见证点 ----
    pi = arb.pi()
    w = {}
    for q in (1, 3):
        lo, _ = point_stats(q, arb(0), arb(0), C)
        _, hi = point_stats(q, pi, pi, C)
        w[q] = (lo[2], hi[3])
    assert arb("0.160639") < w[1][0] and w[1][1] < arb("6.550813")
    assert arb("0.902848") < w[3][0] and w[3][1] < arb("7.053971")
    print("B 荷一窗端：sup h_(3) in [", w[1][0].str(7), ",", win13[1][0].str(6), "]，inf h_(4) in [",
          win13[1][1].str(6), ",", w[1][1].str(7), "]")
    print("B 荷三窗端：sup h_(3) in [", w[3][0].str(7), ",", win13[3][0].str(6), "]，inf h_(4) in [",
          win13[3][1].str(6), ",", w[3][1].str(7), "]")

    overlap = {2: ((3763, 1748), (1102, 1098)), 4: ((582, 2771), (3059, 3075)),
               6: ((700, 3113), (921, 887))}
    expect = {2: ("16.259", "14.575"), 4: ("16.281", "14.621"), 6: ("16.153", "14.690")}
    for q, (p1, p2) in overlap.items():
        lo, _ = point_stats(q, grid(p1[0], 4096), grid(p1[1], 4096), C)
        _, hi = point_stats(q, grid(p2[0], 4096), grid(p2[1], 4096), C)
        assert lo[2] > hi[3]
        assert lo[2] > arb(expect[q][0]) and hi[3] < arb(expect[q][1])
        print("B 荷%d：h_(3)(p1) > %s > %s > h_(4)(p2)，p1=2pi%s/4096，p2=2pi%s/4096"
              % (q, expect[q][0], expect[q][1], p1, p2))

    # 联合相位点 (u,v)=(e^{iX},e^{iY})，X,Y in 2pi Z/64；荷 q 读 (U,V)=(u^q,v^q)。
    fail_pts = {1: (29, 35), 2: (1, 31), 3: (33, 26), 4: (37, 37), 5: (0, 0), 6: (15, 14),
                7: (60, 0), 8: (42, 43), 9: (6, 10), 10: (2, 34), 11: (54, 22), 12: (42, 21),
                13: (41, 23), 14: (42, 21), 16: (10, 42), 17: (54, 42), 18: (48, 48)}
    cache = {}

    def at64(q, ij):
        if (q, ij) not in cache:
            cache[q, ij] = point_stats(q, grid(q*ij[0], 64), grid(q*ij[1], 64), C)
        return cache[q, ij]
    margins = {}
    for r, ij in fail_pts.items():
        # 存在荷 q1、q2 使 h_{q1,(r)} 的下界大于 h_{q2,(r+1)} 的上界（逐项比较精确端点）。
        best = None
        for q1 in CHARGES:
            for q2 in CHARGES:
                d = at64(q1, ij)[0][r-1]-at64(q2, ij)[1][r]
                if d > 0 and (best is None or d > best):
                    best = d
        assert best is not None, r
        margins[r] = best
    assert arb("0.3799") < margins[14] < arb("0.3800")
    assert all(margins[r] > arb("1.383") for r in margins if r != 14)
    print("B 联合点见证余量:", {r: m.str(4) for r, m in margins.items()})
    print("B 秩 r in {1..18}\\{15}：各在一个联合点上 max_q h_(r) > min_q h_(r+1)+0.3799，不是全环面共同分界")
    lo00 = {q: at64(q, (0, 0)) for q in CHARGES}
    assert all(lo00[q][1][2] < arb("7.53559") for q in CHARGES)
    assert all(lo00[q][0][3] > arb("10.6225") for q in CHARGES)
    print("B U=V=1 处五荷共同三模分界成立：max_q h_(3) < 7.53559 < 10.6225 < min_q h_(4)")

    # ---- C. 单步容许秩 {1,3,7,15} 下的狄拉克海窗 ----
    A, B = arb(13), arb(17)
    cells = {q: adaptive(q, A, B, C) for q in CHARGES}
    lo6, hi6 = point_stats(6, grid(3456, 4096), grid(1033, 4096), C)
    lo1, hi1 = point_stats(1, arb(0), arb(0), C)
    wa23 = (lo6[1]+lo6[2])/2
    wa47 = (hi1[3]+hi1[4]+hi1[5]+hi1[6])/4
    assert wa23 > arb("11.68872") and wa47 < arb("17.97989")
    assert sum(cells.values()) == 32608
    assert arb("17.97989")-arb("11.68872") < arb("6.2912")
    print("C 格心数", cells, "；max a23 in [", wa23.str(7), ", 13)，min a47 in (17,", wa47.str(7), "]")

    # 原生能级（delta 取遍 [0,delta_max]）在容许秩 {1,3,7,15} 下的选择。
    scales, _ = prev["native_scales"](arb(DMAX/2, DMAX/2))
    a23_6 = (hi6[1]+hi6[2])/2       # 用于在见证点上的完整比较
    a47_6 = (lo6[3]+lo6[4]+lo6[5]+lo6[6])/4
    a815_6 = sum(lo6[7:15])/8
    a23_1 = (hi1[1]+hi1[2])/2
    a815_1 = sum(lo1[7:15])/8

    def choose(mu):
        if mu < A and mu > B:
            raise AssertionError
        if mu > A and mu < B:
            return 3
        if mu < wa23:
            # 荷六见证点：a23>mu，又 a47,a815>mu，故 E(1) 唯一最小。
            assert a47_6 > mu and a815_6 > mu
            return "1@q6"
        if mu > wa47:
            # 荷一 U=V=1：a23<mu<a815 且 a47<mu：E(7)-E(3)<0，E(15)-E(7)>0，E(7)-E(1)<0。
            assert a23_1 < mu and a815_1 > mu
            assert 2*(a23_1-mu)+4*(wa47-mu) < 0
            return "7@q1"
        return None
    chosen = {k: choose(v) for k, v in scales.items()}
    for k, v in chosen.items():
        if k.startswith(("站点度 (10", "站点方差 (10", "最小割", "最小度", "固定")):
            assert v == "1@q6", (k, v)
        elif k.startswith(("站点度 (18", "站点方差 (18", "最大")):
            assert v == "7@q1", (k, v)
        else:
            assert v == 3, (k, v)
    print("C 原生能级：")
    for k, v in scales.items():
        print("   ", k, v.str(8), "->", chosen[k])

    # ---- D. 前 100 枚实际正零点 ----
    heights = acb.zeta_zeros(1, 100)
    lap = shell["laplacian"]
    no_cut, t_fail, ok18 = 0, [], 0
    for n, z in enumerate(heights, start=1):
        g = z.imag
        delta = arb(23)/25*arb(2)**-20/(g**2+shell["soft"]["point"]["study"]["offset"])
        cz = 15-delta/3
        T = 9+(1-delta)**2
        st = {}
        for q in CHARGES:
            ev = [v.real for v in lap(q, g).eig(multiple=True)]
            st[q] = order_stats([hrange(x.lower(), x.upper(), cz) for x in ev])
        # 五荷共同三模分界不存在：存在 q1,q2 使 h_{q1,(3)} 下界 > h_{q2,(4)} 上界。
        assert any(st[q1][0][2] > st[q2][1][3] for q1 in CHARGES for q2 in CHARGES), n
        no_cut += 1
        # {1,3,7,15} 下 mu=18：逐荷 a23 上界 < 18 < a47 下界。
        assert all((st[q][1][1]+st[q][1][2])/2 < 18 for q in CHARGES), n
        assert all(sum(st[q][0][3:7])/4 > 18 for q in CHARGES), n
        ok18 += 1
        # mu=T：若某荷 a23 下界 > T 则在该零点不选 3；否则逐荷 a23 上界 < T 且 a47 下界 > T。
        if any((st[q][0][1]+st[q][0][2])/2 > T for q in CHARGES):
            t_fail.append(n)
        else:
            assert all((st[q][1][1]+st[q][1][2])/2 < T for q in CHARGES), n
            assert all(sum(st[q][0][3:7])/4 > T for q in CHARGES), n
    assert no_cut == 100 and ok18 == 100
    assert len(t_fail) == 67 and t_fail[:5] == [1, 3, 4, 5, 6] and 2 not in t_fail
    print("D 前 100 枚零点：无一有五荷共同两步三模分界；{1,3,7,15} 下 mu=18 处处选 3；mu=T 失败于",
          len(t_fail), "枚：", t_fail)
    print("用时 %.1f 秒" % (time.time()-t0))


if __name__ == "__main__":
    verify()
