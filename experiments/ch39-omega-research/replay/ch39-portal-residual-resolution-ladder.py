"""三次谱门户留出残差的分辨率阶梯：T 的包围半径来源、随步数的缩放，以及两枚未判定候选的判定。

配套说明：门户留出残差的四个未判定项在1048576步内全部判定而按一阶区间Euler的实测半径律证明可行还需多十四个以上数量级的步数.md。

T_n 的 Arb 球由两部分构成：
  (a) 规范流的区间 Euler 积分（带先验盒），半径约与步数成反比；
  (b) 交点反推逆耦合 1/alpha_s 的包围，经 T 对逆耦合的真实变化传入。
(b) 用交点像映射 G 的不动点迭代 X <- X ∩ G(X) 收紧：盒内交点 x=G(x) 必在 G(X) 中，交集仍包住它，
交点分支约定同 ch39-gauge-window-candidate-scan.py。

运行方式：
  python3 <本脚本>                  只核验数据文件中各球之间的关系与全部断言，并现场重算 22027 的一级抽查（约 7 秒）；
  python3 <本脚本> --check N        重算零点 N 的数据（不写文件），逐球与数据文件比对：相同或相交才通过
                                    （拟合与未判定的六枚每枚约 3 分钟，如 12610 约 192 秒；已判定八枚每枚约 18 秒）；
  python3 <本脚本> --generate N     重算六枚之一并写入数据文件（覆盖该枚条目）；
  python3 <本脚本> --holdout N      重算八枚已判定候选之一在 131072 步、原逆耦合盒下的数据并写入数据文件。
数据文件为同目录的 ch39-portal-residual-resolution-ladder-data.json；Arb 球以精确的二进中心与半径存储。
"""
import json
import sys
import time
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).parent
DATA = ROOT / "ch39-portal-residual-resolution-ladder-data.json"
prev = run_path(str(ROOT / "ch39-rotation-character-portal.py"))
hold = run_path(str(ROOT / "ch39-cubic-portal-holdout.py"))
common, chiral, fixed, two = prev["common"], prev["chiral"], prev["fixed"], prev["two"]
sm = fixed["sm"]
SCALE, UNIT, EPS = prev["SCALE"], prev["UNIT"], prev["TOLERANCE"]

FIT = (11661, 11934, 12207, 28325)
OPEN = (12610, 22027)
LADDER = FIT + OPEN
DECIDED = tuple(h for h in hold["HOLD"] if h not in OPEN)
assert set(DECIDED) == set(hold["CUBIC_DECIDED"])
FLOW_STEPS = (65536, 131072, 262144, 524288, 1048576)
CROSS_STEPS = (32768, 65536, 131072)
CROSS_ITER = 3


def enc(a):
    (mm, me), (rm, re) = a.mid().man_exp(), a.rad().man_exp()
    return [str(mm), int(me), str(rm), int(re)]


def dec(v):
    mid = arb(int(v[0]))*arb(2)**v[1]
    return mid+arb(0, arb(int(v[2]))*arb(2)**v[3])


def levels_of(model, read, n):
    rows = read(acb.zeta_zero(n).imag)
    return {q: chiral["ordered"]([h for h, z in rows if z == q]) for q in (0, *model["charges"])}


def table_of(model, levels, start):
    return [(start+l, w, B, f, s) for l, w, B, f, s in two["thresholds"](model, levels, False)]


def image(model, initial, table, L, inverse, steps):
    """复写 crossing_enclosure 的像映射（无额外质量、不重定中心）。"""
    pi = arb.pi()
    definite = [r for r in table if r[0] < L]
    possible = [r for r in table if not (r[0] > L)]
    segments = fixed["gauge_segments"](L, possible, steps)
    corr, _ = two["flow"]([*initial[:2], inverse], segments, sm["D"], sm["C"], min_inverse=8)
    rhs, mat = [], []
    for i in range(3):
        b0 = sm["BETA"][i]+sum(r[1][i] for r in definite)
        a = initial[i] if i != 2 else 0
        val = 2*pi*(a+corr[i])+sum(r[1][i]*r[0] for r in definite)
        val -= sum(w[i]*(L-l).max(0) for l, w, _, _, _ in possible[len(definite):])
        rhs.append([val])
        mat.append([b0, 2*pi, -2*pi*int(i == 2)])
    return arb_mat(mat).solve(arb_mat(rhs)).entries()


def generate(n):
    model, _, _, initial, _ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    S, bases = prev["rotation"]()
    # candidate 断言规范窗口判据并给出 W、x；其 T 用原逆耦合盒与 65536 步。
    alpha, T0, W, x, _, _, _ = prev["candidate"](model, initial, read, S, bases, n)
    ctx.prec = 192
    entry = {"W": [enc(w) for w in W], "x": [enc(v) for v in x], "T_candidate_65536": enc(T0)}
    levels = levels_of(model, read, n)
    if n in DECIDED:
        ctx.prec = 192
        T, _ = common["portal_at"](model, initial, (levels, 1/alpha, table_of(model, levels, (UNIT/model["mz"]).log())),
                                   SCALE, steps=131072)
        entry["T_original_131072"] = enc(T)
        return entry
    ctx.prec = 192
    mu, _, alpha_c = fixed["crossing_enclosure"](model, list(initial), levels, False, 724, "29.040", ".1184")
    ctx.prec = 192
    start = (arb("2.435e18")/724**2/model["mz"]).log()
    ctable = table_of(model, levels, start)
    L, inverse = (mu/model["mz"]).log(), 1/alpha_c
    entry["inverse_original"] = enc(inverse)
    cross = []
    for steps in CROSS_STEPS:
        for _ in range(CROSS_ITER):
            r = image(model, initial, ctable, L, inverse, steps)
            assert r[0].overlaps(L) and r[2].overlaps(inverse)
            L, inverse = L.intersection(r[0]), inverse.intersection(r[2])
        cross.append({"steps": steps, "inverse": enc(inverse)})
        print(n, "交点步数", steps, "逆耦合半径", arb(inverse.rad()).str(4), flush=True)
    entry["cross"] = cross
    ptable = table_of(model, levels, (UNIT/model["mz"]).log())
    flow = []
    for steps in FLOW_STEPS:
        ctx.prec = 192
        t = time.time()
        T, _ = common["portal_at"](model, initial, (levels, inverse, ptable), SCALE, steps=steps)
        flow.append({"steps": steps, "T": enc(T)})
        print(n, "规范流步数", steps, T.str(15), "半径", arb(T.rad()).str(4), round(time.time()-t, 1), flush=True)
    entry["flow"] = flow
    if n == 22027:
        # 半径拆分：中心逆耦合点上的流半径；工作精度的影响。
        point = arb(inverse.mid())
        entry["point"] = []
        for steps in FLOW_STEPS[:4]:
            ctx.prec = 192
            T, _ = common["portal_at"](model, initial, (levels, point, ptable), SCALE, steps=steps)
            entry["point"].append({"steps": steps, "T": enc(T)})
        entry["precision"] = []
        for prec in (128, 192, 256):
            ctx.prec = prec
            T, _ = common["portal_at"](model, initial, (levels, point, ptable), SCALE, steps=FLOW_STEPS[0])
            entry["precision"].append({"prec": prec, "T": enc(T)})
        ctx.prec = 192
        entry["point_inverse"] = enc(point)
    return entry


def same_balls(old, new, path=""):
    """逐球比对：四元组视为 Arb 球，要求相交；其余字段要求相同。"""
    if isinstance(old, list) and len(old) == 4 and isinstance(old[1], int) and isinstance(old[0], str):
        a, b = dec(old), dec(new)
        assert a.overlaps(b), path
        return 1
    if isinstance(old, dict):
        assert set(old) == set(new), path
        return sum(same_balls(old[k], new[k], f"{path}/{k}") for k in old)
    if isinstance(old, list):
        assert len(old) == len(new), path
        return sum(same_balls(a, b, f"{path}[{i}]") for i, (a, b) in enumerate(zip(old, new)))
    assert old == new, path
    return 0


def store(n, entry):
    data = json.loads(DATA.read_text()) if DATA.exists() else {}
    data[str(n)] = entry
    DATA.write_text(json.dumps(data, indent=1, sort_keys=True))


def holdout(rows, targets, row, target):
    k = len(rows)
    D = arb_mat(rows).det()
    assert not D.contains(0)
    lam = [arb_mat([row if j == i else rows[j] for j in range(k)]).det()/D for i in range(k)]
    weight = sum(abs(v) for v in lam)
    return target-sum(v*t for v, t in zip(lam, targets))+arb(0, (EPS*weight).upper()), weight, lam


def verify():
    ctx.prec = 192
    raw = json.loads(DATA.read_text())
    assert set(map(int, raw)) == set(LADDER) | set(DECIDED)
    W = {int(k): [dec(v) for v in e["W"]] for k, e in raw.items()}
    x = {int(k): [dec(v) for v in e["x"]] for k, e in raw.items()}

    # 一、交点逆耦合的收紧：每级包含于上一级与原始盒，半径随交点步数近似反比下降。
    inv_final = {}
    for n in LADDER:
        e = raw[str(n)]
        prior = dec(e["inverse_original"])
        rads = [arb(prior.rad())]
        for c in e["cross"]:
            box = dec(c["inverse"])
            assert prior.contains(box)
            prior = box
            rads.append(arb(box.rad()))
        assert rads[1]/rads[-1] > arb(3)  # 32768 到 131072 步至少收紧三倍
        assert rads[0]/rads[-1] > arb(10)  # 比原始盒至少紧十倍
        assert rads[-1] < arb("2.6e-5")
        # 交点步数每加倍，逆耦合半径缩小 1.9 至 2.1 倍（与步数成反比）。
        assert all(arb("1.9") < u/v < arb("2.1") for u, v in zip(rads[1:], rads[2:]))
        inv_final[n] = prior
        print("逆耦合包围半径：原始、32768、65536、131072 步", n, [r.str(3) for r in rads])

    # 二、T 的阶梯：相邻两级嵌套；半径乘步数近似常数（流部分与步数成反比）。
    T = {n: [dec(f["T"]) for f in raw[str(n)]["flow"]] for n in LADDER}
    for n in LADDER:
        assert [f["steps"] for f in raw[str(n)]["flow"]] == list(FLOW_STEPS)
        for a, b in zip(T[n], T[n][1:]):
            assert a.contains(b)
        # 65536 步下，收紧后的逆耦合盒给出的 T 含于原 candidate 的 T 之内。
        assert dec(raw[str(n)]["T_candidate_65536"]).contains(T[n][0])
        products = [arb(t.rad())*s for t, s in zip(T[n], FLOW_STEPS)]
        assert all(arb("4.4e24") < p < arb("6.4e24") for p in products), [p.str(3) for p in products]
        print("T 半径 × 步数", n, [p.str(3) for p in products])

    # 半径拆分（22027）：逆耦合取中心点时流半径 × 步数在 (4.40,4.47)e24 内（2^16 至 2^19 步）；精度 128–256 位不改变半径。
    e = raw["22027"]
    point = [arb(dec(p["T"]).rad())*p["steps"] for p in e["point"]]
    assert all(arb("4.40e24") < p < arb("4.47e24") for p in point)
    box_minus_point = [arb(t.rad())-arb(dec(p["T"]).rad()) for t, p in zip(T[22027], e["point"])]
    assert all(arb("1.6e18") < d < arb("1.85e18") for d in box_minus_point)
    prec = [arb(dec(p["T"]).rad()) for p in e["precision"]]
    assert all(abs(p/prec[1]-1) < arb("1e-3") for p in prec)
    print("22027 中心点流半径×步数", [p.str(4) for p in point], "逆耦合盒附加半径", [d.str(3) for d in box_minus_point],
          "精度 128/192/256 位的半径", [p.str(4) for p in prec])

    # 三、各级的留出残差。
    rows = lambda ns, k: [x[n][:k] for n in ns]
    results = {}
    for level, steps in enumerate(FLOW_STEPS):
        tf = [T[n][level] for n in FIT]
        out = {}
        for h in OPEN:
            r, weight, lam = holdout(rows(FIT, 4), tf, x[h][:4], T[h][level])
            out[("cubic", h)] = (r, weight)
        tri = (11934, 12207, 28325)
        for name, k, src in (("quad", 3, x), ("char", 3, W)):
            r, weight, _ = holdout([src[n][:k] for n in tri], [T[n][level] for n in tri], src[12610][:k],
                                   T[12610][level])
            out[(name, 12610)] = (r, weight)
        results[steps] = out
        print("步数", steps, {f"{a}-{b}": (r.lower().str(4), r.upper().str(4)) for (a, b), (r, _) in out.items()})

    # 各级残差包围嵌套（同一组 lambda，T 嵌套）。
    for a, b in zip(FLOW_STEPS, FLOW_STEPS[1:]):
        for key in results[a]:
            assert results[a][key][0].contains(results[b][key][0])

    # 各项首个不含零的步数、符号与最细一级的下界。
    expect = {("cubic", 22027): (-1, 524288, "4.4e19"), ("cubic", 12610): (1, 262144, "4.3e19"),
              ("quad", 12610): (1, 131072, "1.1e20"), ("char", 12610): (-1, 1048576, "4.4e18")}
    for key, (sign, first, low) in expect.items():
        for steps in FLOW_STEPS:
            r = results[steps][key][0]
            if steps < first:
                assert r.contains(0), (key, steps)
            else:
                assert sign*r > 0, (key, steps)
        assert sign*results[FLOW_STEPS[-1]][key][0] > arb(low)
        lo_hi = {("cubic", 22027): ("-1.066e20", "-4.42e19"), ("cubic", 12610): ("4.36e19", "6.90e19"),
                 ("quad", 12610): ("1.10e20", "1.355e20"), ("char", 12610): ("-2.942e19", "-4.49e18")}[key]
        final = results[FLOW_STEPS[-1]][key][0]
        assert arb(lo_hi[0]) < final.lower() and final.upper() < arb(lo_hi[1])
        print("判定", key, "首个不含零的步数", first, "最细一级", results[FLOW_STEPS[-1]][key][0].str(5))

    # 12610 的 |r| 上界序列单调下降；最细一级的包围在 (4.36,6.90)e19 内、严格为正。
    ups = [abs(results[s][("cubic", 12610)][0]).upper() for s in FLOW_STEPS]
    assert all(a > b for a, b in zip(ups, ups[1:]))
    bands = (("2.02e20", "2.03e20"), ("1.31e20", "1.32e20"), ("9.55e19", "9.56e19"), ("7.77e19", "7.78e19"),
             ("6.89e19", "6.90e19"))
    assert all(arb(a) < u < arb(b) for u, (a, b) in zip(ups, bands))
    print("12610 三次残差 |r| 上界序列", [arb(u).str(4) for u in ups])
    # 十枚留出候选中最小的 |r| 下界属于 12610，约为 eps 的 5e15 倍以上。
    low12610 = abs(results[FLOW_STEPS[-1]][("cubic", 12610)][0]).lower()
    assert arb("4.3e19") < low12610 < abs(results[FLOW_STEPS[-1]][("cubic", 22027)][0]).lower()
    assert all(low12610 < arb(v) for _, v in hold["CUBIC_DECIDED"].values())
    assert low12610/EPS > arb("5e15")
    # 逆耦合盒附加半径与流半径相等的步数：在 2^21 与 2^22 之间（实测系数外推）。
    assert arb("4.40e24")/arb("1.85e18") > 2**21 and arb("4.47e24")/arb("1.6e18") < 2**22

    # 四、已判定八枚：用收紧后的拟合 T（含于原包围）重算，符号与原下界不变。
    for h in DECIDED:
        Th = dec(raw[str(h)]["T_original_131072"])
        r, _, _ = holdout(rows(FIT, 4), [T[n][1] for n in FIT], x[h][:4], Th)
        sign, low = hold["CUBIC_DECIDED"][h]
        assert sign*r > arb(low)
        print("已判定留出候选", h, "残差", r.str(5))

    # 五、可行性一侧。可行指存在 |e_i|<=eps 使 |r_h(e)|<=eps，等价于 |c|<=eps(1+w)，
    # c=T_h-sum lambda_i T_i，w=sum|lambda_i|。要由区间证书证明它，须 rad(c)<=eps(1+w)；
    # 而 rad(c)>=(1+w)·min_n rad(T_n)，故每枚 T 的半径须小于 eps，与 w 无关。
    # 按一阶区间 Euler 的实测流半径系数（22027 盒心点、2^16 至 2^19 步，>4.40e24）外推所需步数。
    need = arb("4.40e24")/EPS
    assert arb("5.6e20") < need < arb("5.7e20")
    assert need/FLOW_STEPS[-1] > arb("5.3e14")
    print("证明可行时每枚 T 的流半径须小于 eps，所需步数（一阶区间 Euler 的半径律外推）>", need.lower().str(3),
          "，是最细一级的", (need/FLOW_STEPS[-1]).lower().str(3), "倍")

    # 现场抽查：22027 在 65536 步、收紧逆耦合盒下重算 T，与数据文件中的球相交。
    model, _, _, initial, _ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    levels = levels_of(model, read, 22027)
    ptable = table_of(model, levels, (UNIT/model["mz"]).log())
    live, _ = common["portal_at"](model, initial, (levels, inv_final[22027], ptable), SCALE, steps=FLOW_STEPS[0])
    assert live.overlaps(T[22027][0]) and abs(arb(live.rad())/arb(T[22027][0].rad())-1) < arb("1e-6")
    print("现场重算 22027 的 65536 步 T", live.str(12))


if __name__ == "__main__":
    if len(sys.argv) == 3 and sys.argv[1] in ("--generate", "--holdout"):
        n = int(sys.argv[2])
        assert (n in LADDER) if sys.argv[1] == "--generate" else (n in DECIDED)
        store(n, generate(n))
    elif len(sys.argv) == 3 and sys.argv[1] == "--check":
        n = int(sys.argv[2])
        assert n in LADDER or n in DECIDED
        old = json.loads(DATA.read_text())[str(n)]
        count = same_balls(old, json.loads(json.dumps(generate(n))))
        print("重算与数据文件逐球比对通过", n, "球数", count)
    else:
        verify()
