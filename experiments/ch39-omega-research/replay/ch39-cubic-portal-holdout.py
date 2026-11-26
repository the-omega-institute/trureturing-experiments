"""四枚规范候选唯一确定的三次谱门户、以及三枚确定的二次谱门户和轮换群代数门户，在留出候选上的残差。

配套说明：四枚规范候选唯一确定的三次谱门户在十枚留出候选中八枚上的残差被证书认定非零.md。
T 的 Arb 球是对连续规范流的严格包围（区间 Euler 格式加先验盒，离散误差计入半径），
步数只影响半径、不改变被包围的量；脚本断言 131072 步的球含于 65536 步的球。

读数与 ch39-rotation-character-portal.py 相同：一圈对角门户读数加纯电弱两圈有限阈值 T_n，
抵消条件 sum_k a_k x_k(n) = T_n（谱多项式）或 sum_chi p_chi W_chi(n) = T_n（群代数），容差 eps=125.25^2/2 GeV^2。
k 维族在 k 枚拟合零点上的解唯一时，留出零点 h 上的残差对 T 线性：
    r_h = T_h - sum_i lambda_i T_i,  lambda_i = det(把第 i 行换成 h 的行) / det(拟合矩阵)。
拟合零点上只要求 |残差| <= eps 时，r_h 再加宽 eps*sum|lambda_i|。
k+1 枚零点上的族内极小极大残差为 |y.T|/||y||_1（y 为带号子式）；当 k 枚之外只多一枚 h 时它等于 |r_h|/(1+sum|lambda_i|)。
"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

ROOT = Path(__file__).parent
prev = run_path(str(ROOT / "ch39-rotation-character-portal.py"))
scan = run_path(str(ROOT / "ch39-gauge-window-candidate-scan.py"))
common, chiral = prev["common"], prev["chiral"]

FIT = prev["ZEROS"]
assert FIT == (11661, 11934, 12207, 28325)
# 第 12208 至 13700 枚的全部规范窗口候选（逐编号排查见 ch39-gauge-window-candidate-scan.py），
# 以及研究档案中出现过、高度大于 10000 的另外五枚（档案处方不同，这里由 candidate 按同一判据重新认证；
# 未按规则穷尽档案中的全部此类零点）。
SCANNED = scan["HITS"]
NAMED = (19872, 22027, 48067, 58214, 95478)
HOLD = SCANNED + NAMED
STATE = {"steps": 131072}
COARSE = 65536
EPS = prev["TOLERANCE"]

# 原 candidate 以 65536 步跑 T 的规范流；这里只把步数换成 STEPS，其余判据与读数不变。
_portal_at = common["portal_at"]


def _finer(*args, steps=None, **kwargs):
    return _portal_at(*args, steps=STATE["steps"], **kwargs)


common["portal_at"] = _finer

# 活跃十四层的特征标序列。A 为 11661 的序列，B 为 11934 的序列。
CLASS_A = [1, 2, 0, 2, 0, 1, 0, 0, 2, 1, 2, 1, 1, 0]
CLASS_B = [1, 2, 0, 2, 1, 0, 0, 2, 1, 0, 2, 1, 1, 0]
CLASS_C = [1, 2, 0, 2, 1, 0, 0, 0, 2, 1, 2, 1, 1, 0]
CLASS_D = [1, 2, 0, 2, 1, 0, 0, 1, 2, 0, 2, 1, 1, 0]
CLASSES = {11661: "A", 12207: "A", 28325: "A", 12336: "A", 12886: "A", 13699: "A", 58214: "A",
           11934: "B", 12610: "B", 48067: "B", 95478: "B", 19872: "C", 22027: "C", 13291: "D"}

# 三次门户的留出残差：已判定者给出符号与 |r| 的下界，未判定者给出 |r| 的上界（GeV^2）。
CUBIC_DECIDED = {12336: (1, "3.6e21"), 12886: (1, "3.3e21"), 13291: (1, "8.1e21"), 13699: (1, "4.2e21"),
                 19872: (-1, "1.4e20"), 48067: (1, "4.1e20"), 58214: (1, "8.5e21"), 95478: (1, "6.6e19")}
CUBIC_OPEN = {12610: "1.9e20", 22027: "3.7e20"}


def gather():
    model, _, _, initial, _ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    S, bases = prev["rotation"]()
    data = {}
    for n in FIT + HOLD:
        # candidate 断言：五荷三模、0.117<alpha_s<0.119、Q 以下恰有十四个活跃层、显式算子与 S 对易。
        alpha, T, W, x, pattern, _, _ = prev["candidate"](model, initial, read, S, bases, n)
        ctx.prec = 192
        assert pattern == {"A": CLASS_A, "B": CLASS_B, "C": CLASS_C, "D": CLASS_D}[CLASSES[n]]
        assert T > 0 and T.rad() < arb("7e19") and all(w > 0 for w in W)
        # 矩与特征标权的包围远比 T 紧：残差包围几乎全部来自 T。
        assert all(v.rad() < arb("1e-40")*v for v in [*W, *x])
        if n in (19872, 95478):
            # 这两枚的判定依赖 131072 步；核对两种步数的包围嵌套。
            STATE["steps"] = COARSE
            _, coarse, _, _, _, _, _ = prev["candidate"](model, initial, read, S, bases, n)
            STATE["steps"] = 131072
            ctx.prec = 192
            assert coarse.contains(T) and coarse.rad() > T.rad()
            print("步数嵌套", n, "65536 步", coarse.str(8), "含 131072 步", T.str(8), flush=True)
        data[n] = (T, W, x)
        print("零点、强耦合、T/GeV²、T 包围半径、特征标序列类", n, alpha.str(6), T.str(8),
              arb(T.rad()).str(3), CLASSES[n], flush=True)
    return data


def holdout(rows, targets, row, target):
    k = len(rows)
    D = arb_mat(rows).det()
    assert not D.contains(0)
    lam = [arb_mat([row if j == i else rows[j] for j in range(k)]).det()/D for i in range(k)]
    weight = sum(abs(v) for v in lam)
    r = target-sum(v*t for v, t in zip(lam, targets))+arb(0, (EPS*weight).upper())
    return r, weight, lam


def minimax(rows, targets):
    k = len(rows)
    y = [(-1)**i*arb_mat([r for j, r in enumerate(rows) if j != i]).det() for i in range(k)]
    return abs(sum(a*t for a, t in zip(y, targets)))/sum(abs(a) for a in y)


def verify():
    data = gather()
    T = {n: data[n][0] for n in data}
    span_lo, span_hi = T[FIT[0]], T[FIT[0]]
    for n in data:
        span_lo, span_hi = span_lo.min(T[n]), span_hi.max(T[n])
    assert arb("1.0451e26") < span_lo and span_hi < arb("1.0475e26")
    assert span_hi-span_lo > arb("2.2e23")
    assert EPS/span_lo < arb("7.6e-23")
    print("十四枚候选 T 的范围/GeV²", span_lo.str(6), span_hi.str(6))

    # 平凡基准：只拟合一枚的常数门户 a0=T_f/x0_f 在其余各枚上的相对残差。
    base, worst_by_fit = arb(0), {}
    for f in data:
        m = arb(0)
        for h in data:
            if h != f:
                m = m.max(abs(T[h]-T[f]/data[f][2][0]*data[h][2][0])/T[h])
        worst_by_fit[f] = m
        base = base.max(m)
    assert arb("8.4e-3") < base < arb("8.5e-3")
    # 拟合零点选得最好时（11661），在其余十三枚上仍有相对残差 >4.9e-3。
    assert arb("4.9e-3") < worst_by_fit[11661] < arb("5.0e-3")
    assert all(worst_by_fit[11661] < v for f, v in worst_by_fit.items() if f != 11661)
    base_best = worst_by_fit[11661]
    print("单枚常数门户：任一拟合零点、任一其他零点的相对残差上界", base.upper().str(4),
          "最佳拟合零点 11661 的最坏相对残差", worst_by_fit[11661].str(5))

    # 其余留出候选离四枚拟合零点都远。
    assert abs(T[12610]-T[11934]) < arb("1.8e21")
    for h in HOLD:
        if h != 12610:
            assert all(abs(T[h]-T[f]) > arb("1.4e22") for f in FIT)
    print("|T_12610-T_11934|/GeV² 上界", abs(T[12610]-T[11934]).upper().str(3))

    # 一、四枚定出的三次谱门户在十枚留出候选上的残差。
    rows = [data[n][2][:4] for n in FIT]
    targets = [T[n] for n in FIT]
    worst_rel = arb(0)
    least_five = None
    lows = {}
    for h in HOLD:
        r, weight, lam = holdout(rows, targets, data[h][2][:4], T[h])
        if h == 12610:
            # 12610 的矩行几乎等于 11934 的矩行。
            assert abs(lam[1]-1) < arb(".02") and sum(abs(lam[i]) for i in (0, 2, 3)) < arb(".1")
            assert all(abs(lam[i]) < arb(".05") for i in (0, 2, 3))
            print("12610 的 lambda", [v.str(4) for v in lam])
        five = abs(r)/(1+weight)
        rel = abs(r)/T[h]
        worst_rel = worst_rel.max(rel)
        if h in CUBIC_DECIDED:
            sign, low = CUBIC_DECIDED[h]
            assert sign*r > arb(low)
            assert (abs(r) / EPS) > arb("8e15")
            assert five > 0
            least_five = five if least_five is None else least_five.min(five)
            lows[h] = abs(r).lower()
            print("三次门户留出残差/GeV²（已判定）", h, r.str(5), "五枚极小极大下界", five.lower().str(4))
        else:
            assert r.contains(0) and abs(r) < arb(CUBIC_OPEN[h])
            print("三次门户留出残差/GeV²（包围含零）", h, r.str(5))
    assert set(CUBIC_DECIDED) | set(CUBIC_OPEN) == set(HOLD)
    assert worst_rel < arb("9e-5")
    # 三次门户的最坏相对残差夹在已判定者的下界与全体上界之间；
    # 与最佳单枚常数门户相比，它只小五十五至六十二倍。
    worst_low = arb(0)
    for h in CUBIC_DECIDED:
        worst_low = worst_low.max(arb(lows[h])/T[h])
    assert arb(55) < base_best.lower()/worst_rel.upper() and base_best.upper()/worst_low.lower() < arb(62)
    assert least_five > arb("1.36e19")
    # 八枚已判定者中，|r| 的下界最小的是 95478。
    assert all(lows[95478] < v for h, v in lows.items() if h != 95478)
    print("三次门户留出残差相对 T 的上界", worst_rel.upper().str(4),
          "已判定八枚的五枚极小极大残差下界", least_five.lower().str(4))

    # 二、三枚定出的二次谱门户与群代数门户，在其余十一枚上的残差。
    for name, row, bound in (("二次谱门户", lambda n: data[n][2][:3], "1.65e20"),
                             ("群代数门户", lambda n: data[n][1], "1.15e20")):
        undecided, least = [], None
        for triple in combinations(FIT, 3):
            for h in HOLD+tuple(m for m in FIT if m not in triple):
                r, _, _ = holdout([row(n) for n in triple], [T[n] for n in triple], row(h), T[h])
                if r.contains(0):
                    undecided.append((triple, h))
                    assert abs(r) < arb("3e20")
                else:
                    least = abs(r) if least is None else least.min(abs(r))
        assert undecided == [((11934, 12207, 28325), 12610)], undecided
        assert least > arb(bound)
        print(name, "四组三枚拟合 × 十一枚留出：43 项已判定，最小 |残差| 下界/GeV²", least.lower().str(4),
              "唯一未判定（不计为失败）", undecided)

    # 三、族层面：十四枚中任取 k+1 枚的族内极小极大残差。
    names = sorted(data)
    for name, row, size, open_count, bound in (
            ("群代数门户", lambda n: data[n][1], 4, 26, "6.7e18"),
            ("二次谱门户", lambda n: data[n][2][:3], 4, 52, "3.27e18"),
            ("三次谱门户", lambda n: data[n][2][:4], 5, 279, "4.9e17")):
        undecided, least, total = 0, None, 0
        for chosen in combinations(names, size):
            total += 1
            m = minimax([row(n) for n in chosen], [T[n] for n in chosen])
            if m.lower() > 0:
                least = m if least is None else least.min(m)
            else:
                undecided += 1
        assert undecided == open_count and least > arb(bound)
        assert arb(bound)/EPS > arb("6e13")
        print(name, "十四枚中", total, "个", size, "元子集：未判定（不计为失败）", undecided, "，其余排除，最小残差下界/GeV²",
              least.lower().str(4))


if __name__ == "__main__":
    verify()
