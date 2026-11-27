"""十四枚规范候选上同时抵消两圈希格斯阈值的谱多项式门户的最低次数。

配套说明：十四枚规范候选上谱多项式门户的最低次数在6109个子集上被证书认定恰为候选数减一而全体十四枚只证得次数在7与13之间.md。

读数与 ch39-rotation-character-portal.py、ch39-cubic-portal-holdout.py 相同：
    x_k(n) = sum_i w_i z_i^k，z=(m/80M)^2，w=m^2(1+2log(Q/m))/(16 pi^2)，对 Q 以下十四个荷三活跃层求和；
    次数 d 的谱多项式门户 (a_0..a_d) 在候选集 S 上可行，指 max_{n in S} |sum_k a_k x_k(n) - T_n| <= eps。
d(S) 为可行的最低次数。对 |S|=N：
    R_d(S) = min_a max_n |x(n).a - T_n|；N x (N-1) 矩阵列满秩时 R_{N-2}(S) = |y.T|/||y||_1，y 为带号 (N-1) 阶子式；
    N x N 矩矩阵（次数 0..N-1）非奇异时次数 N-1 的门户唯一且精确抵消，故 d(S) <= N-1；
    R_{N-2}(S) > eps 时 d(S) >= N-1（次数更低的族是其子族）。
T 取已入库的最紧包围：拟合四枚与 12610、22027 取 ch39-portal-residual-resolution-ladder-data.json 中
2^20 步、收紧逆耦合盒的一级；其余八枚取同一文件中 131072 步、原逆耦合盒的一级。

运行：python3 <本脚本>            重算并与数据文件中的被排除子集逐一比对
      python3 <本脚本> --write    重算并写出数据文件
"""
import json
import sys
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).parent
OUT = ROOT / "ch39-spectral-portal-degree-count-data.json"
lad = run_path(str(ROOT / "ch39-portal-residual-resolution-ladder.py"))
prev, hold = lad["prev"], lad["hold"]
soft = prev["soft"]
UNIT, SCALE, EPS = prev["UNIT"], prev["SCALE"], prev["TOLERANCE"]
CLASS = {"A": hold["CLASS_A"], "B": hold["CLASS_B"], "C": hold["CLASS_C"], "D": hold["CLASS_D"]}
NAMES = tuple(sorted(hold["FIT"] + hold["HOLD"]))
assert len(NAMES) == 14
PREC = 320
DEG = 14


def layers(n, S, bases):
    """复写 candidate 中显式荷三算子的谱部分：活跃层的特征标与 g=m/M。"""
    ctx.prec = PREC
    height = acb.zeta_zero(n).imag
    unit_matrix = prev["charge_matrix"](3, height)
    base, neutral, _ = soft["vector_pencil"](height)
    amount, star = prev["spoke_removal"](S, height)
    mat = unit_matrix-amount*star
    trace = (base+neutral).trace().real
    assert trace == 270
    center = (trace-6*amount)/18
    labelled = []
    for chi, B in enumerate(bases):
        for eig in (B.conjugate().transpose()*mat*B).eig(multiple=True, algorithm="rump"):
            assert eig.imag.contains(0)
            v = eig.real
            labelled.append((v*abs(v-center), chi))
    assert len(labelled) == 19
    labelled.sort(key=lambda t: float(t[0].mid()))
    assert all(a[0] < b[0] for a, b in zip(labelled, labelled[1:]))
    active = []
    for g, chi in labelled[3:]:
        mass = UNIT*g
        assert mass < SCALE or mass > SCALE
        if mass < SCALE:
            active.append((g, chi))
    assert len(active) == 14
    return active


def gather():
    ctx.prec = PREC
    S, bases = prev["rotation"]()
    raw = json.loads(lad["DATA"].read_text())
    dec = lad["dec"]
    data = {}
    for n in NAMES:
        active = layers(n, S, bases)
        ctx.prec = PREC
        assert [c for _, c in active] == CLASS[hold["CLASSES"][n]]
        x, W = [arb(0)]*DEG, [arb(0)]*3
        for g, chi in active:
            mass = UNIT*g
            w = mass**2*(1+2*(SCALE/mass).log())/(16*arb.pi()**2)
            W[chi] += w
            for k in range(DEG):
                x[k] += w*(g/80)**(2*k)
        e = raw[str(n)]
        # 与 candidate 写入数据文件的前四阶矩、三特征标权同一。
        assert all(a.overlaps(dec(v)) for a, v in zip(x, e["x"]))
        assert all(a.overlaps(dec(v)) for a, v in zip(W, e["W"]))
        assert all(v.rad() < arb("1e-45")*v for v in [*x, *W])
        top = dec(e["T_candidate_65536"])
        if "flow" in e:
            levels = [dec(f["T"]) for f in e["flow"]]
            assert [f["steps"] for f in e["flow"]] == list(lad["FLOW_STEPS"])
            assert top.contains(levels[0]) and all(a.contains(b) for a, b in zip(levels, levels[1:]))
            T = levels[-1]
            assert T.rad() < arb("6.1e18")
        else:
            T = dec(e["T_original_131072"])
            assert top.contains(T) and T.rad() < arb("6.4e19")
        assert T > 0
        data[n] = (T, x, W)
    return data


def holder(rows, targets):
    """N x (N-1) 情形的精确极小极大值 |y.T|/||y||_1 的包围；同时返回 y。"""
    N = len(rows)
    y = [(-1)**i*arb_mat([r for j, r in enumerate(rows) if j != i]).det() for i in range(N)]
    l1 = sum((abs(a) for a in y), arb(0))
    return abs(sum(a*t for a, t in zip(y, targets)))/l1, y


def mask(S):
    return sum(1 << NAMES.index(n) for n in S)


def verify(write):
    data = gather()
    T = {n: data[n][0] for n in NAMES}
    x = {n: data[n][1] for n in NAMES}
    for n in NAMES:
        assert x[n][0] > 0 and T[n] > arb("1.04e26")
    # 最小的阈值半径属于六枚阶梯候选，超过 eps 的 7e14 倍：任何 R 的包围都无法认定 R <= eps。
    radmin = min(arb(T[n].rad()) for n in NAMES)
    assert all(arb(T[n].rad()) > arb("5.9e18") for n in NAMES)
    assert all(arb(T[n].rad()) > arb("5e19") for n in NAMES if n not in lad["LADDER"])
    assert radmin/EPS > arb("7e14")
    worst_y = arb(0)

    certified = {}
    excluded = {N: [] for N in range(1, 15)}
    least = {}
    for N in range(1, 15):
        for S in combinations(NAMES, N):
            # 定理 1：N x N 矩矩阵非奇异。
            D = arb_mat([x[n][:N] for n in S]).det()
            assert not D.contains(0)
            if N == 1:
                # 零门户（“次数 -1”）的残差是 T 本身。
                ok = True
                value = T[S[0]]
            else:
                value, y = holder([x[n][:N-1] for n in S], [T[n] for n in S])
                # 列满秩：y 的分量是 (N-1) 元子集的矩行列式，已在上一层认证非零。
                assert all(not a.contains(0) for a in y)
                for a in y:
                    worst_y = worst_y.max(arb(a.rad())/abs(a))
                # 相容性核对：包围的上端总超过 eps。
                assert value.upper() > EPS
                # 盒内取值论证：T_n 取遍 [c_n-r_n, c_n+r_n] 时 sup|y.T| = |y.c| + sum|y_n| r_n，
                # 它除以 ||y||_1 仍超过 eps，故只依赖现有包围的证书推不出 R_{N-2}(S) <= eps。
                cy = abs(sum(a*arb(T[n].mid()) for a, n in zip(y, S)))
                ry = sum((abs(a)*arb(T[n].rad()) for a, n in zip(y, S)), arb(0))
                l1 = sum((abs(a) for a in y), arb(0))
                assert (cy+ry)/l1 > EPS and ry/l1 > arb(radmin.lower())*arb("0.999999")
                ok = value.lower() > EPS
            certified[mask(S)] = ok
            if ok:
                excluded[N].append(mask(S))
                least[N] = value if N not in least else least[N].min(value)
            if N == 14:
                norms = arb(1)
                for n in S:
                    norms *= sum((v**2 for v in x[n]), arb(0)).sqrt()
                hadamard = abs(D)/norms
    counts = [len(excluded[N]) for N in range(1, 15)]
    total = [len(list(combinations(NAMES, N))) for N in range(1, 15)]
    print("各 N 的子集数", total)
    print("各 N 上被证书认定 d(S)=N-1 的子集数", counts)
    assert counts == [14, 91, 364, 991, 1875, 1581, 1162, 31, 0, 0, 0, 0, 0, 0]
    # 措辞的比例：四至六元层过半，七元层约三分之一。
    assert all(2*counts[N-1] > total[N-1] for N in (4, 5, 6))
    assert 0.33 < counts[6]/total[6] < 0.34
    assert sum(counts) == 6109 and sum(total) == 16383
    assert [t-c for t, c in zip(total, counts)] == [0, 0, 0, 10, 127, 1422, 2270, 2972, 2002, 1001, 364, 91, 14, 1]
    assert sum(total[8:]) == 3473 and sum(total)-sum(counts) == 10274
    assert worst_y < arb("1e-22")
    print("子式 y 的最坏相对半径", worst_y.upper().str(3))
    bounds = {2: "2.2e21", 3: "3.2e19", 4: "7.1e17", 5: "6.0e17", 6: "8.6e16", 7: "4.3e16", 8: "2.7e17"}
    for N, b in bounds.items():
        assert least[N] > arb(b) and arb(b)/EPS > arb("5e12")
        print("N =", N, "被排除子集中 R_{N-2} 的最小下界/GeV²", least[N].lower().str(4))
    # 比较的是各层最小下界（球的下端），不是球本身。
    assert all(arb(least[7].lower()) < arb(least[N].lower()) for N in bounds if N != 7)
    # 全体十四枚矩矩阵的 Hadamard 比（|det| 除以各行二范数之积）：非奇异但极近退化。
    assert arb("4.2e-115") < hadamard < arb("4.3e-115")
    print("十四枚 14x14 矩矩阵的 Hadamard 比", hadamard.str(4))

    # 未判定子集的显式名单（N=4、5 两层）。
    undecided4 = [S for S in combinations(NAMES, 4) if not certified[mask(S)]]
    assert len(undecided4) == 10
    for S in undecided4:
        print("N=4 未判定（不计为失败）", S)
    assert undecided4 == [(11661, 13699, 19872, 28325), (11661, 19872, 22027, 48067), (11934, 12610, 13291, 95478),
                          (11934, 12610, 19872, 22027), (11934, 12610, 48067, 95478), (11934, 13699, 48067, 58214),
                          (12207, 12886, 13291, 95478), (12610, 13699, 48067, 58214), (12886, 13699, 48067, 95478),
                          (13291, 19872, 48067, 95478)]
    # 三十一个被认定的八元子集都含 11661、12207、13291、58214。
    eight = [{n for i, n in enumerate(NAMES) if m >> i & 1} for m in excluded[8]]
    common8 = set.intersection(*eight)
    assert common8 == {11661, 12207, 13291, 58214}
    print("八元被认定子集的公共成员", sorted(common8))

    # 由单调性得到的下界 L(S)：d(S) >= max{|S'|-1 : S' ⊆ S 已认定}。
    full = (1 << 14)-1
    L = {}
    for m in range(1, full+1):
        size = bin(m).count("1")
        if certified[m]:
            L[m] = size-1
        else:
            L[m] = max(L[m & ~(1 << i)] for i in range(14) if m >> i & 1)
    assert L[full] == 7
    gaps = {}
    for m in range(1, full+1):
        size = bin(m).count("1")
        gaps.setdefault(size, set()).add(size-1-L[m])
    assert {N: sorted(v) for N, v in gaps.items()} == {
        1: [0], 2: [0], 3: [0], 4: [0, 1], 5: [0, 1], 6: [0, 1, 2], 7: [0, 1, 2, 3], 8: [0, 1, 2, 3],
        9: [1, 2, 3, 4], 10: [2, 3, 4, 5], 11: [3, 4, 5], 12: [4, 5], 13: [5, 6], 14: [6]}
    print("全体十四枚：7 <= d <= 13；各 N 的上下界差集合", {N: sorted(v) for N, v in gaps.items()})

    # 群代数门户与已有子集计数：收紧 T 后重数四元子集。
    ga = sum(1 for S in combinations(NAMES, 4)
             if holder([data[n][2] for n in S], [T[n] for n in S])[0].lower() > EPS)
    assert ga == 987
    print("群代数门户四元子集被排除数", ga, "/ 1001")

    # 仅数值：以包围中心代替 T 时全体十四枚的 R_12，及其包围半径的加权平均。
    rows = [x[n][:13] for n in NAMES]
    _, y = holder(rows, [T[n] for n in NAMES])
    l1 = sum((abs(a) for a in y), arb(0))
    mid = abs(sum(a*arb(T[n].mid()) for a, n in zip(y, NAMES)))/l1
    res = sum((abs(a)*arb(T[n].rad()) for a, n in zip(y, NAMES)), arb(0))/l1
    assert arb("3.3e17") < mid < arb("3.5e17") and arb("3.1e19") < res < arb("3.2e19")
    assert arb(92) < res/mid < arb(93)
    assert EPS == arb("7843.78125")
    print("仅数值：十四枚 R_12 的中心读数/GeV²", mid.str(4), "加权半径/GeV²", res.str(4))

    record = {"names": list(NAMES), "excluded_masks": {str(N): excluded[N] for N in range(1, 15)}}
    if write:
        OUT.write_text(json.dumps(record, indent=0)+"\n")
        print("已写出", OUT.name)
    else:
        stored = json.loads(OUT.read_text())
        assert stored == record
        print("与数据文件逐一一致：", OUT.name)


if __name__ == "__main__":
    verify("--write" in sys.argv[1:])
