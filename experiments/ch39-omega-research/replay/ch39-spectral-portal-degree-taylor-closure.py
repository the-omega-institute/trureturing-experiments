"""十四枚规范候选上谱多项式门户最低次数的完全判定：高阶区间 Taylor 规范流与显式谱把阈值包围收紧到容差以下。

配套说明：十四枚规范候选的全部16383个子集上谱多项式门户的最低次数都被证书认定恰为候选数减一而全体十四枚在次数不超过十七的单项式支撑中至少要十三项.md。

读数、门户族、容差与候选同 ch39-spectral-portal-degree-count.py。阈值 T_n 的包围在两处收紧：
  (a) 两圈规范流（三个逆耦合修正与三个 Yukawa 对数，六维常微分方程）改用 ORDER 阶区间 Taylor 法：
      每步先用 Picard 判据 u_k + [0,h] f([t_k,t_k+h], B) ⊆ B 证明解留在先验盒 B 内，
      再取 u(t_k+h) ∈ sum_{j<p} u_j(t_k,u_k) h^j + h^p u_p([t_k,t_k+h], B)（Lagrange 余项），
      Taylor 系数由 1/v 与 exp 的幂级数递推精确给出；接口与 ch39-two-tier-dirac-running.py 的
      integrate_segments 相同，替换后交点像映射与 portal_at 不作其他改动；
  (b) 阈值谱不用 spectral_supplier 的 1e-12 宽带，改用显式算子 charge_matrix(q) - a*star
      （三条辐边对每个荷都是单位传输，ch39-unit-seventy-five-crossing.py 的 shell_spectra 注释）
      按轮换特征标分块的 Rump 本征值包围；脚本断言每个显式球含于原宽带球中。
交点包围从 crossing_enclosure 的同一种子盒出发，再迭代 X <- X ∩ G(X)；交点分支约定同前。

运行：python3 <本脚本>        全部重算（单核约四分钟）
"""
import sys
import time
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).parent
dc = run_path(str(ROOT / "ch39-spectral-portal-degree-count.py"))
lad, prev = dc["lad"], dc["prev"]
common, chiral, fixed, two = prev["common"], prev["chiral"], prev["fixed"], prev["two"]
soft, matter = prev["soft"], prev["matter"]
NAMES, EPS, UNIT, SCALE = dc["NAMES"], dc["EPS"], dc["UNIT"], dc["SCALE"]
ORDER, HMAX = 18, arb("0.25")
SPARSE_DEGREE = 17


# ---------- (a) 区间 Taylor 规范流 ----------
def coeffs(a, s, c0, l0, tau0, M, D, C, p):
    """在 (tau0, c0, l0) 处解的 0..p 阶 Taylor 系数；逆耦合 v_i = a_i - s_i tau/(2pi) + c_i，g=1/v，y=exp(2l)。"""
    pi = arb.pi()
    c, l = [[v] for v in c0], [[v] for v in l0]
    v = [[a[i]-s[i]*tau0/(2*pi)+c0[i]] for i in range(3)]
    g = [[1/v[i][0]] for i in range(3)]
    w = [[2*l0[i]] for i in range(3)]
    y = [[w[i][0].exp()] for i in range(3)]
    k8, k32, k16 = 1/(8*pi**2), 1/(32*pi**3), 1/(16*pi**2)
    for k in range(p):
        gk, yk = [g[j][k] for j in range(3)], [y[j][k] for j in range(3)]
        # 右端对 g、y 线性：与 integrate_segments 的 dx、dy 逐项相同。
        dc_ = [-sum(M[i][j]*gk[j] for j in range(3))*k8+sum(D[i][j]*yk[j] for j in range(3))*k32 for i in range(3)]
        tr = 3*(yk[0]+yk[1])+yk[2]
        diff = (yk[0]-yk[1], yk[1]-yk[0], yk[2])
        dl = [(arb(3)/2*diff[i]+tr-4*pi*sum(C[i][j]*gk[j] for j in range(3)))*k16 for i in range(3)]
        m = k+1
        for i in range(3):
            c[i].append(dc_[i]/m)
            l[i].append(dl[i]/m)
            v[i].append(c[i][m]+(-s[i]/(2*pi) if m == 1 else 0))
            g[i].append(-sum(v[i][j]*g[i][m-j] for j in range(1, m+1))/v[i][0])
            w[i].append(2*l[i][m])
            y[i].append(sum(j*w[i][j]*y[i][m-j] for j in range(1, m+1))/m)
    return c, l, g


def horner(cs, h):
    r = cs[-1]
    for t in reversed(cs[:-1]):
        r = r*h+t
    return r


def taylor_flow(initial, segments, D, C, record=None, *, min_inverse=8):
    assert record is None
    pi, p = arb.pi(), ORDER
    correction = [arb(0) for _ in range(3)]
    logs = [arb(x).log() for x in (".98", ".016", ".010")]
    baseline = initial[:]
    for span, coeff, matrix, _ in segments:
        n = max(1, int((span.upper()/HMAX).ceil().unique_fmpz()))
        h = span/n
        for k in range(n):
            tau0 = k*h
            Bt = tau0+arb(0).union(1)*h
            cs, ls, _ = coeffs(baseline, coeff, correction, logs, tau0, matrix, D, C, 1)
            Bc = [correction[i]+arb(0).union(1)*h*cs[i][1] for i in range(3)]
            Bl = [logs[i]+arb(0).union(1)*h*ls[i][1] for i in range(3)]
            ok = False
            for _ in range(30):
                Bc = [b+arb(0, arb(b.rad())*arb("0.2")+arb("1e-40")) for b in Bc]
                Bl = [b+arb(0, arb(b.rad())*arb("0.2")+arb("1e-40")) for b in Bl]
                cb, lb, gb = coeffs(baseline, coeff, Bc, Bl, Bt, matrix, D, C, 1)
                assert min_inverse > 0 and all(1/gb[i][0] > min_inverse for i in range(3))
                ec = [correction[i]+arb(0).union(1)*h*cb[i][1] for i in range(3)]
                el = [logs[i]+arb(0).union(1)*h*lb[i][1] for i in range(3)]
                # Picard 判据：成立时解在 [t_k, t_k+h] 上存在且留在 B 内。
                if all(Bc[i].contains(ec[i]) and Bl[i].contains(el[i]) for i in range(3)):
                    ok = True
                    break
                Bc = [b.union(e) for b, e in zip(Bc, ec)]
                Bl = [b.union(e) for b, e in zip(Bl, el)]
            assert ok
            cb, lb, _ = coeffs(baseline, coeff, Bc, Bl, Bt, matrix, D, C, p)
            cp, lp, _ = coeffs(baseline, coeff, correction, logs, tau0, matrix, D, C, p-1)
            hp = h**p
            correction = [horner(cp[i], h)+hp*cb[i][p] for i in range(3)]
            logs = [horner(lp[i], h)+hp*lb[i][p] for i in range(3)]
        baseline = [a-s*span/(2*pi) for a, s in zip(baseline, coeff)]
    return correction, [v+c for v, c in zip(baseline, correction)]


# ---------- (b) 显式谱 ----------
def explicit_levels(model, read, S, bases, n):
    ctx.prec = 192
    height = acb.zeta_zero(n).imag
    old = lad["levels_of"](model, read, n)
    amount, star = prev["spoke_removal"](S, height)
    # 辐边减量的 Weyl 范数界 4a 远小于 spectral_supplier 所留的 1e-12 宽带。
    assert height > 11000 and 4*amount < arb("2.8e-14")
    base, neutral, _ = soft["vector_pencil"](height)

    def spec(M):
        out = []
        for B in bases:
            for e in (B.conjugate().transpose()*M*B).eig(multiple=True, algorithm="rump"):
                assert e.imag.contains(0)
                out.append(e.real)
        assert len(out) == 19
        return out
    lap = base+neutral
    assert all(z.contains(0) for z in (prev["charge_matrix"](0, height)-lap).entries())
    neu = spec(lap-amount*star)
    assert sum(1 for v in neu if v.contains(0)) == 1
    neu = [v for v in neu if not v.contains(0)]
    assert len(neu) == 18 and all(v > arb(".01") for v in neu)
    center = sum(neu)/18
    rows = [(v*abs(v-center), 0) for v in neu]
    for q in model["charges"]:
        M = prev["charge_matrix"](q, height)
        assert all(z.contains(0) for z in (M-base-matter["nonspokes"](q, height)).entries())
        rows += [(v*abs(v-center), q) for v in spec(M-amount*star)]
    new = {q: chiral["ordered"]([h for h, z in rows if z == q]) for q in (0, *model["charges"])}
    # 每个显式顺序统计量含于 spectral_supplier 的 1e-12 宽带之中。
    for q in new:
        assert len(new[q]) == len(old[q]) and all(a.contains(b) for a, b in zip(old[q], new[q]))
    assert all(arb(h.rad()) < arb("1e-45") for q in new for h in new[q])
    return new


def threshold(model, initial, levels):
    ctx.prec = 192
    mu, _, alpha_c = fixed["crossing_enclosure"](model, list(initial), levels, False, 724, "29.040", ".1184")
    ctx.prec = 192
    ctable = lad["table_of"](model, levels, (arb("2.435e18")/724**2/model["mz"]).log())
    L, inverse = (mu/model["mz"]).log(), 1/alpha_c
    rads = [arb(inverse.rad())]
    for it in range(40):
        r = lad["image"](model, initial, ctable, L, inverse, 0)
        assert r[0].overlaps(L) and r[2].overlaps(inverse)
        L, inverse = L.intersection(r[0]), inverse.intersection(r[2])
        rads.append(arb(inverse.rad()))
        if it > 3 and not rads[-1] < rads[-2]*arb("0.9"):
            break
    ptable = lad["table_of"](model, levels, (UNIT/model["mz"]).log())
    T, _ = common["portal_at"](model, initial, (levels, inverse, ptable), SCALE, steps=0)
    return T, inverse


def main():
    global ORDER, HMAX
    t0 = time.time()
    data = dc["gather"]()  # 旧包围、四阶以内的矩与 data 文件的一致性
    model, _, _, initial, _ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 192
    read = chiral["source"]["spectral_supplier"](model)
    S, bases = prev["rotation"]()
    two["flow"] = taylor_flow
    T, LEV = {}, {}
    for n in NAMES:
        LEV[n] = explicit_levels(model, read, S, bases, n)
        T[n], inv = threshold(model, initial, LEV[n])
        old = data[n][0]
        # 新包围含于已入库的最紧包围（一阶 Euler，2^20 或 131072 步）。
        assert old.contains(T[n]) and T[n] > arb("1.04e26")
        assert arb("5.9e18") < arb(old.rad()) < arb("6.4e19")
        assert arb(T[n].rad()) < arb("3e-4") and arb(old.rad())/arb(T[n].rad()) > arb("2e22")
        print(n, "T =", T[n].str(22), "半径", arb(T[n].rad()).str(3), "原半径/新半径",
              (arb(old.rad())/arb(T[n].rad())).str(3), "逆耦合半径", arb(inv.rad()).str(3), flush=True)
    radmax = max(arb(T[n].rad()) for n in NAMES)
    assert all(arb(T[n].rad()) < arb("3e-4") for n in NAMES) and EPS/arb("3e-4") > arb("2.6e7")
    print("十四枚 T 的半径全部小于 3e-4 GeV²，不到 eps 的 1/(2.6e7)；最大", radmax.str(4), "用时", round(time.time()-t0))

    # 积分格式自洽：换阶数与步长重算一枚，各球两两相交；半径随阶数与步长按 Taylor 余项变化。
    alt = {}
    for order, step in ((10, "0.125"), (14, "0.25"), (14, "0.125")):
        ORDER, HMAX = order, arb(step)
        alt[(order, step)], _ = threshold(model, initial, LEV[22027])
        print("22027 以", order, "阶、步长", step, "重算", alt[(order, step)].str(22), "半径",
              arb(alt[(order, step)].rad()).str(3))
    ORDER, HMAX = 18, arb("0.25")
    balls = [T[22027], *alt.values()]
    assert all(a.overlaps(b) for a in balls for b in balls)
    assert arb(alt[(10, "0.125")].rad()) > arb(1e4) and arb(alt[(14, "0.25")].rad()) > arb(100)
    assert arb(alt[(14, "0.125")].rad()) < arb(1)
    # 对照：同一积分器若沿用 spectral_supplier 的 1e-12 宽带谱，T 的半径留在远大于 eps 的地板上。
    Tw, _ = threshold(model, initial, lad["levels_of"](model, read, 22027))
    assert Tw.contains(T[22027]) and arb(Tw.rad()) > arb("1e14")
    print("22027 以宽带谱重算的 T 半径", arb(Tw.rad()).str(3))

    # ---------- 全部子集的最低次数 ----------
    ctx.prec = 320
    act = {n: dc["layers"](n, S, bases) for n in NAMES}
    ctx.prec = 320
    X = {}
    for n in NAMES:
        x = [arb(0)]*(SPARSE_DEGREE+1)
        for g, chi in act[n]:
            mass = UNIT*g
            w = mass**2*(1+2*(SCALE/mass).log())/(16*arb.pi()**2)
            for k in range(SPARSE_DEGREE+1):
                x[k] += w*(g/80)**(2*k)
        assert all(a.overlaps(b) for a, b in zip(x, data[n][1]))
        X[n] = x
    holder = dc["holder"]
    counts = [0]*15
    least = {}
    for N in range(1, 15):
        for Sub in combinations(NAMES, N):
            assert not arb_mat([X[n][:N] for n in Sub]).det().contains(0)
            if N == 1:
                v = T[Sub[0]]
            else:
                v, y = holder([X[n][:N-1] for n in Sub], [T[n] for n in Sub])
                assert all(not a.contains(0) for a in y)
            if v.lower() > EPS:
                counts[N] += 1
                least[N] = v if N not in least else least[N].min(v)
    total = [len(list(combinations(NAMES, N))) for N in range(1, 15)]
    print("各 N 上被认定 d(S)=N-1 的子集数", counts[1:])
    assert counts[1:] == total and sum(counts) == 16383
    lows = {N: least[N].lower() for N in range(2, 15)}
    for N in range(2, 15):
        print("N =", N, "R_{N-2} 最小下界/GeV²", arb(lows[N]).str(4), "与 eps 之比", (arb(lows[N])/EPS).str(3))
    table = {2: "2.2e21", 3: "7.1e19", 4: "1.38e17", 5: "3.43e17", 6: "2.54e16", 7: "4.89e15", 8: "2.17e15",
             9: "4.39e13", 10: "7.18e13", 11: "6.47e14", 12: "8.93e15", 13: "1.04e15", 14: "2.13e15"}
    assert all(arb(lows[N]) > arb(b) and arb(lows[N]) < 2*arb(b) for N, b in table.items())
    floor = min(arb(v) for v in lows.values())
    assert arb("4.39e13") < floor < arb("4.40e13") and floor/EPS > arb("5.6e9")
    floorN = [N for N in range(2, 15) if not arb(lows[N]) > arb("1.5e14")]
    assert floorN == [9, 10] and arb(lows[9]) < arb(lows[10])
    assert least[9].upper() < arb("4.40e13") and least[9].upper() < least[10].lower()
    print("最小下界所在层", floorN, arb(floor).str(4))
    # 分辨率门槛：R_{N-2}(S)=|y.T|/||y||_1 对 T 的上确界范数 1-Lipschitz，故任何半径都小于 (L_N-eps)/2、
    # 且包住真值的阈值包围都能认定第 N 层全部子集（充分条件）。
    gate = {N: (arb(lows[N])-EPS)/2 for N in range(2, 15)}
    print("各层门槛 (L_N-eps)/2 / GeV²", {N: g.str(3) for N, g in gate.items()})
    assert arb("2.19e13") < gate[9] < arb("2.20e13") and all(gate[N] >= gate[9] for N in gate)
    assert arb("6.9e16") < gate[4] < arb("7.0e16") and arb("1.2e16") < gate[6] < arb("1.3e16")
    # 计数说明所用的旧半径（均大于 5.9e18，见上文断言）高于 N>=4 各层的门槛。
    assert all(gate[N] < arb("5.9e18") for N in range(4, 15))
    full, _ = holder([X[n][:13] for n in NAMES], [T[n] for n in NAMES])
    assert arb("2.1395e15") < full < arb("2.1396e15")
    # 计数说明实测读数 6：以旧包围中心代替 T 的中心读数，与认证值相差两个数量级以上。
    _, yfull = holder([X[n][:13] for n in NAMES], [T[n] for n in NAMES])
    l1 = sum((abs(a) for a in yfull), arb(0))
    centre = abs(sum(a*arb(data[n][0].mid()) for a, n in zip(yfull, NAMES)))/l1
    assert arb("3.3e17") < centre < arb("3.5e17") and arb("3.3e17")/full.upper() > 100
    print("旧包围中心读数", centre.str(4), "与认证值之比下界", (arb("3.3e17")/full.upper()).str(4))
    # 群代数门户（三个特征标权）在 1001 个四元子集上用新 T 重算。
    gres = [holder([data[n][2] for n in Sub], [T[n] for n in Sub])[0] for Sub in combinations(NAMES, 4)]
    ga = sum(1 for v in gres if v.lower() > EPS)
    gmin = min(arb(v.lower()) for v in gres)
    assert ga == 1001 and gmin > arb("4.3e18")
    print("群代数门户四元子集被排除数", ga, "/ 1001，最小下界", gmin.str(4))
    print("全体十四枚 R_12 ∈", full.str(10))

    # ---------- 全体十四枚的稀疏支撑 ----------
    # 若 m<=12 项支撑 K（指数取自 0..17）在全体上可行，则任一含 K 的 12 项支撑可行，
    # 且在每个十三元子集上可行；十三元子集上 R_K = |y.T|/||y||_1（y 为 12 阶子式且全非零时）。
    ts = time.time()
    supports = list(combinations(range(SPARSE_DEGREE+1), 12))
    worst = None
    for K in supports:
        hit = False
        for drop in NAMES:
            Sp = [n for n in NAMES if n != drop]
            v, y = holder([[X[n][k] for k in K] for n in Sp], [T[n] for n in Sp])
            if any(a.contains(0) for a in y):
                continue
            if v.lower() > EPS:
                hit = True
                worst = v if worst is None else worst.min(v)
                break
        assert hit, K
    assert len(supports) == 18564 and worst.lower() > arb("5.5e12")
    print("次数不超过 17 的 18564 个十二项支撑全部被排除；各支撑所用见证的最小下界/GeV²", arb(worst.lower()).str(4),
          "用时", round(time.time()-ts))
    print("总用时", round(time.time()-t0))


if __name__ == "__main__":
    main()
