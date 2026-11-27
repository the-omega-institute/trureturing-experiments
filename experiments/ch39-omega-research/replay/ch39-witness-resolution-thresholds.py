"""四问模型内见证的分辨率门槛：零点高度球与两道钟相位球上的区间证书，以及盒外的有理失败点。"""
from collections import defaultdict
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from runpy import run_path
import time

from flint import acb, acb_mat, arb, ctx, fmpz_mat

ROOT = Path(__file__).resolve().parent
phys = run_path(str(ROOT/"ch39-first-three-mode-observables.py"))
uni = run_path(str(ROOT/"ch39-unit-seventy-five-crossing.py"))
smooth, soft = uni["smooth"], uni["soft"]
matter = smooth["matter"]
geo, pole = matter["geo"], matter["pole"]
joint = phys["joint"]
chiral, vacuum, gauge = joint["chiral"], joint["vacuum"], joint["gauge"]
two = gauge["two"]
model, initial = phys["model_inputs"]()
CHARGES = model["charges"]
OFFSET = soft["point"]["study"]["offset"]
ctx.prec = 128
# 两道钟相位为 A_d·γ，A_d=θ(d)-θ(244)。
A = {d: geo["theta"](d)-geo["theta"](244) for d in (724, 804)}
CLOCK = [(i, j, d, z) for i, j, d, z in geo["edges"] if d in (724, 804) and pole not in (i, j)]
_base, _neutral, _ = soft["vector_pencil"](arb(10001))
BASE = _base
_int = fmpz_mat([[int(t.real.unique_fmpz()) for t in row] for row in (_base+_neutral).tolist()])
NEUTRAL = [z.real for z, k in _int.charpoly().complex_roots() if z != 0 for _ in range(k)]
assert len(NEUTRAL) == 18


# ---------- 带荷谱：任意给定的钟相位，加上相位误差球 ----------
def nonspokes(q, ph):
    m = acb_mat(19, 19)
    for i, j, d, z in geo["edges"]:
        if pole in (i, j):
            continue
        link = (z*ph.get(d, 1))**q
        m[i, i] += 1
        m[j, j] += 1
        m[i, j] = -link
        m[j, i] = -link.conjugate()
    return m


def weyl(q, delta):
    """相位误差 |ε_d|≤δ_d 时扰动矩阵的谱范数上界（厄米矩阵取最大绝对行和）。"""
    rows = [arb(0)]*19
    for i, j, d, _ in CLOCK:
        e = (abs(q)*delta[d]).min(arb(2))
        rows[i] += e
        rows[j] += e
    r = arb(0)
    for x in rows:
        r = r.max(x)
    return r.upper()


def eigen_enclosure(q, ph, delta):
    """对相位误差球内每个厄米矩阵 H 包住其排序本征值。
    Weyl：|λ_i(H)-λ_i(H0)|≤W。孤立本征值再用固定试向量 v 的 Temple–Kato 不等式
    ρ-η²/(b-ρ)≤λ_i≤ρ+η²/(ρ-a)，ρ=v*Hv/v*v，η=‖Hv-ρv‖/‖v‖，(a,b) 由 Weyl 区间给出。"""
    H0 = BASE+nonspokes(q, ph)
    ev, Rv = H0.eig(right=True)
    lam = chiral["ordered"]([z.real for z in ev])
    order = sorted(range(19), key=lambda i: float(ev[i].real.mid()))
    W = weyl(q, delta)
    E = acb_mat(19, 19)
    for i, j, d, z in CLOCK:
        e = -(z*ph[d])**q*(acb(0, q*arb(0, delta[d])).exp()-1)
        E[i, j] = e
        E[j, i] = e.conjugate()
    H = H0+E
    out = []
    for idx, col in enumerate(order):
        enc = lam[idx]+arb(0, W)
        a = lam[idx-1].upper()+W if idx > 0 else None
        b = lam[idx+1].lower()-W if idx < 18 else None
        if (a is None or enc.lower() > a) and (b is None or enc.upper() < b):
            v = acb_mat([[acb(Rv[r, col].real.mid(), Rv[r, col].imag.mid())] for r in range(19)])
            vh = v.transpose().conjugate()
            nv = (vh*v)[0, 0].real
            Hv = H*v
            rho = (vh*Hv)[0, 0].real/nv
            res = Hv-v*acb(rho)
            eta2 = (res.transpose().conjugate()*res)[0, 0].real/nv
            if (a is None or rho > a) and (b is None or rho < b):
                lo = rho if b is None else rho-eta2/(b-rho)
                hi = rho if a is None else rho+eta2/(rho-a)
                enc = enc.intersection(lo.lower().union(hi.upper()))
        out.append(enc)
    return out


def spectra(ph, delta, hlow):
    # 三条辐边的刚度减量 ≤ 4·(23/25)·2^-20/(γ²+offset)，γ≥hlow。
    shift = arb(0).union(4*arb(23)/25*arb(2)**-20/(hlow**2+OFFSET))
    out = {0: [v-shift for v in NEUTRAL]}
    for q in CHARGES:
        out[q] = [v-shift for v in eigen_enclosure(q, ph, delta)]
    return out


def layers_of(levels):
    center = sum(levels[0])/18
    return {q: chiral["ordered"]([v*abs(v-center) for v in eig]) for q, eig in levels.items()}


def union_layers(cells):
    return {q: [_union([c[q][j] for c in cells]) for j in range(len(cells[0][q]))] for q in cells[0]}


def _union(xs):
    u = xs[0]
    for x in xs[1:]:
        u = u.union(x)
    return u


# ---------- 三条见证的判据：True 为全盒认证成立，False 为认证不成立，None 为未判 ----------
def three_mode(L):
    if all(L[q][2] < 9 < L[q][3] for q in CHARGES):
        return True
    if any(L[q][2] > 9 or L[q][3] < 9 for q in CHARGES):
        return False
    return None


def balanced(L):
    """与原配平相同的权与逐质量事件检验，去掉正负质量间隙断言：正权质量向上、负权质量向下取整，
    对球内每个构型累计和都不低于所检验的保守累计和。"""
    rows = vacuum["minimal"]["fields"](L)
    assert all(g % 4 == 0 for g, _, _ in rows if g > 0)
    assert sum(-g for g, _, _ in rows if g < 0) == 958
    assert sum(g//4 for g, _, _ in rows if g > 0) == 162
    residual, events = [], defaultdict(int)
    for g, a, z in rows:
        extra = arb(-g)/6 if g < 0 else arb(g)/12
        weight = int((2*(a+extra)).unique_fmpz())
        residual.append((arb(weight)/2, z))
        if z.is_zero():
            assert weight == 21
            continue
        bound = z*10**15
        cut = bound.upper().ceil() if weight > 0 else bound.lower().floor()
        events[int(cut.unique_fmpz())] += weight
    total, ok = 0, True
    for edge in sorted(events):
        total += events[edge]
        ok = ok and total >= 0
    assert total+21 == 393
    return residual, ok


COEF = (9*arb(3).log()-12*arb(2).log())/2


def signed_mass_gap(L):
    """正权与负权质量包围之间的最小有向间隙：每对取 max(z下-w上, w下-z上)，为负即两球重叠。"""
    rows = vacuum["minimal"]["fields"](L)
    pos = [z for g, _, z in rows if g > 0]
    neg = [z for g, _, z in rows if g < 0]
    out = None
    for z in pos:
        for w in neg:
            v = (z.lower()-w.upper()).max(w.lower()-z.upper())
            out = v if out is None else out.min(v)
    return out


def curvature(L):
    mass = 32*arb("2.435e18")/724**2*L[2][15]
    a = (arb("246.22")**2/(mass*arb("1.22089e19")))**2
    res, ok = balanced(L)
    rows = res+vacuum["low_sector"](4, a)
    bracket = vacuum["newton"](rows, arb(".03")) < 1 < vacuum["newton"](rows, arb(".033"))
    curve = a**2*COEF/(arb.pi()*arb("1e-122"))*(1-arb(0).union(a/(arb(".03")*COEF)))
    return ok and bracket, curve


def curvature_window(L):
    ok, c = curvature(L)
    if ok and arb(".998") < c < arb("1.002"):
        return True
    if c < arb(".998") or c > arb("1.002"):
        return False
    return None


_strict_thresholds = two["thresholds"]


def relaxed_thresholds(model, layers, full):
    """原两圈阈值表，放宽严格次序断言。chiral["ordered"] 把下端点与上端点分别排序，第 j 个球
    包住真实的第 j 个顺序统计量（下断言其两端非降）；同荷结点按下标次序进入表（下断言）；
    不同荷的阈值球严格分开。于是表的次序就是真实排序，逐段 max(0,·) 的包围仍有效。"""
    assert not full
    for q in (1, 3):
        assert all(a.lower() <= b.lower() and a.upper() <= b.upper()
                   for a, b in zip(layers[q], layers[q][1:]))
    knots = []
    for q in (1, 3):
        color = 3 if q == 1 else 1
        indices = [arb(2*color*q*q)/60, arb(color)/2, arb(int(color == 3))]
        casimirs = [arb(q*q)/60, arb(3)/4, arb(4)/3 if color == 3 else 0]
        for j, h in enumerate(layers[q]):
            fermion = int(j >= 3)
            scalar = int(q == 3 and (full or j >= 3))
            if not fermion and not scalar:
                continue
            one = [v*(4*fermion+scalar)/3 for v in indices]
            rep = model["weyl"][model["charges"].index(q)]
            expected = [2*fermion*w+scalar*a for w, a in zip(rep, model["higgs"])]
            assert all((a-b).contains(0) for a, b in zip(one, expected))
            tw = [[indices[i]*(4*(fermion+scalar)*casimirs[k] +
                    ((20*fermion+2*scalar)*arb((0, 2, 3)[i])/3 if i == k else 0))
                   for k in range(3)] for i in range(3)]
            knots.append((h.log(), one, tw, fermion, scalar, q, j))
    knots.sort(key=lambda row: float(row[0].mid()))
    for q in (1, 3):
        order = [k[6] for k in knots if k[5] == q]
        assert order == sorted(order)
    for x in range(len(knots)):
        for y in range(x+1, len(knots)):
            if knots[x][5] != knots[y][5]:
                assert knots[x][0] < knots[y][0]
    return [k[:5] for k in knots]


def alpha(L, relaxed=True):
    two["thresholds"] = relaxed_thresholds if relaxed else _strict_thresholds
    try:
        with redirect_stdout(StringIO()):
            _, _, al = gauge["crossing_enclosure"](model, initial, L, False, 724, "29.04", ".1184", recenter=True)
    finally:
        two["thresholds"] = _strict_thresholds
    return al


def alpha_window(L):
    try:
        al = alpha(L)
    except AssertionError:
        return None
    if arb(".117") < al < arb(".119"):
        return True
    if al < arb(".117") or al > arb(".119"):
        return False
    return None


# ---------- 两种盒 ----------
def height_box(h0, k, m):
    """γ∈[γ_n-2^-k, γ_n+2^-k] 等分为 2^m 格，每格以中心相位加 |A_d|·格半径 的相位球。"""
    r, n = arb(2)**-k, 2**m
    cells = []
    for i in range(n):
        hc = h0+r*(2*i+1-n)/n
        cr = r/n
        ph = {d: acb(0, A[d]*hc).exp() for d in A}
        cells.append(layers_of(spectra(ph, {d: abs(A[d])*cr for d in A}, hc-cr)))
    return union_layers(cells)


def phase_box(h0, x0, y0, k, m):
    """两道钟相位各自取 [ϑ±2^-k]，等分为 2^m×2^m 格；辐边减量按 γ≥γ_n-1 包住。"""
    r, n = arb(2)**-k, 2**m
    cells = []
    for i in range(n):
        for j in range(n):
            ph = {724: acb(0, x0+r*(2*i+1-n)/n).exp(), 804: acb(0, y0+r*(2*j+1-n)/n).exp()}
            cells.append(layers_of(spectra(ph, {724: r/n, 804: r/n}, h0-1)))
    return union_layers(cells)


ZERO = {724: arb(0), 804: arb(0)}


def at_height(g):
    return layers_of(spectra({d: acb(0, A[d]*g).exp() for d in A}, ZERO, g))


def cayley(t):
    """有理 Cayley 参数 t 给出单位圆上的有理点 ((1-t²)+2ti)/(1+t²)。"""
    t = arb(t)
    return acb(1-t*t, 2*t)/(1+t*t)


def cayley_parameter(angle, bits):
    """角 ϑ 附近的二进有理 t≈tan(ϑ/2)，返回 t 及其精确角 2·atan t。"""
    t = (angle/2).tan()
    t = arb(int((t*2**bits).mid().floor().unique_fmpz()))/2**bits
    return t, 2*t.atan()


def at_cayley(h0, tx, ty):
    return layers_of(spectra({724: cayley(tx), 804: cayley(ty)}, ZERO, h0-1))


def dyadic(g, bits):
    return arb(int((g*2**bits).mid().floor().unique_fmpz()))/2**bits


# ---------- 第三问：旗舰零点的味窗只经相对钟角 Δ=ϑ724-ϑ804 ----------
def flavor_tools():
    cp = run_path(str(ROOT/"ch39-three-window-two-probability-cp.py"))["cp"]
    flavor = cp["flavor"]
    data = cp["two"]["spectral_data"]()
    polys = cp["cp_data"](data, flavor["reference_matrix"]())

    def read(c):
        return cp["two"]["probabilities"](c, data)[0][0], cp["cp_value"](c, polys)
    bench = read(flavor["singlet"](11661)[1])
    # 六面平均实部为 a cos Δ + b sin Δ（两钟同转时精确不变，见 singlet_at 的断言）。
    a = flavor["singlet_at"]({724: acb(1), 804: acb(1)})[1]
    b = flavor["singlet_at"]({724: acb(0, 1), 804: acb(1)})[1]
    K, off = (a*a+b*b).sqrt(), acb(a, b).arg()

    def rel(c):
        return [v/w-1 for v, w in zip(read(c), bench, strict=True)]

    def window(c):
        u, j = rel(c)
        if abs(u) < arb("2e-6") and abs(j) < arb("1.2e-6"):
            return True
        if abs(u) > arb("2e-6") or abs(j) > arb("1.2e-6"):
            return False
        return None

    def interval(d0, rho, pieces=64):
        lo, hi = d0-rho, d0+rho
        # 在 (off-π, off) 上 c=K cos(Δ-off) 严格递增，故 c 的值域是两端点之间的区间。
        assert hi-off < 0 and lo-off > -arb.pi()
        clo, chi = K*(lo-off).cos(), K*(hi-off).cos()
        for i in range(pieces):
            cell = (clo+(chi-clo)*arb(i)/pieces).union(clo+(chi-clo)*arb(i+1)/pieces)
            w = window(cell)
            if w is not True:
                return w
        return True
    return flavor, rel, window, interval, K, off


def phase_gap(g1, g2):
    """两枚零点两道钟相位差的最大模（取主值）。"""
    out = arb(0)
    for d in A:
        z = acb(0, A[d]*(g2-g1)).exp()
        out = out.max(abs(z.arg()))
    return out


def log2(x):
    return float((arb(x).log()/arb(2).log()).mid())


def verify_8359(report):
    n = 8359
    ctx.prec = 224
    gamma = acb.zeta_zero(n).imag
    original_radius = gamma.rad()
    ctx.prec = 128
    gamma = acb.zeta_zero(n).imag
    x0, y0 = [acb(0, A[d]*gamma).exp().arg() for d in (724, 804)]
    # 中心点复现原脚本：谱、三模、强耦合、曲率比。
    centre = at_height(gamma)
    ref = uni["shell_spectra"]((0, *CHARGES), zero_index=n)
    levels = spectra({d: acb(0, A[d]*gamma).exp() for d in A}, ZERO, gamma)
    assert all(a.overlaps(b) for q in ref for a, b in zip(levels[q], ref[q]))
    assert three_mode(centre) is True
    ok, curve0 = curvature(centre)
    assert ok and arb(".99935101") < curve0 < arb(".99935103")
    strict, relaxed = alpha(centre, relaxed=False), alpha(centre)
    assert strict.overlaps(relaxed) and arb(".1184112") < relaxed < arb(".1184213")
    assert original_radius < arb(2)**-209
    assert 2**13 <= gamma < 2**14 and arb("8462.069886") < gamma < arb("8462.069887")
    assert arb("0.4121") < A[724] < arb("0.4122") and arb("0.4676") < A[804] < arb("0.4678")
    assert arb("0.055676") < abs(A[724]-A[804]) < arb("0.055677")
    report["8359"] = dict(gamma=gamma, rad=original_radius, alpha=relaxed, curve=curve0)
    print("第 8359 枚：γ", gamma.str(15), "原证书 γ 半径 2^", round(log2(original_radius), 1),
          "强耦合", relaxed.str(10), "曲率比", curve0.str(12))
    return gamma, x0, y0


WITNESSES = (("三模", three_mode), ("曲率", curvature_window), ("强耦合", alpha_window))
# 证书参数：(判据, 高度球 k, 高度格 m, 相位球 k, 相位格 m)
CERT = {"三模": (3, 6, 5, 3), "曲率": (9, 4, 10, 3), "强耦合": (5, 7, 9, 3)}
# 失败点：高度偏移（二进位数）与相位偏移
FAIL_HEIGHT = {"三模": ("0.14", 16), "曲率": ("-0.0034", 24), "强耦合": ("4", 12)}
FAIL_PHASE = {"三模": ("0.0525", "-0.0525"), "曲率": ("-0.0015", "-0.0015"), "强耦合": ("0", "-1.5")}
# 失败点落在 B(k-DROP) 内：DROP=-1 即落在粗一位的盒中，二进门槛恰为 k。
HEIGHT_DROP = {"三模": -1, "曲率": -1, "强耦合": -7}
PHASE_DROP = {"三模": -1, "曲率": -1, "强耦合": -10}
DIST = {"三模": (("0.13999", "0.14"), ("0.05249", "0.05251")),
        "曲率": (("0.0034", "0.0034001"), ("0.00149", "0.00151")),
        "强耦合": (("3.9999", "4"), ("1.49999", "1.50001"))}
NEIGHBOUR = {8358: (("-0.50702", "-0.50701"), ("0.23717", "0.23718"), [False, False, True]),
             8360: (("0.57669", "0.5767"), ("0.26976", "0.26977"), [False, False, True])}


def part_8359(report):
    gamma, x0, y0 = verify_8359(report)
    rows = {}
    for name, test in WITNESSES:
        kh, mh, kp, mp = CERT[name]
        Lh, Lp = height_box(gamma, kh, mh), phase_box(gamma, x0, y0, kp, mp)
        assert test(Lh) is True and test(Lp) is True
        if name == "曲率":
            # 去掉正负质量间隙断言在相位盒上是承重的：两类质量球已重叠，高度盒上仍分开。
            gh, gp = signed_mass_gap(Lh), signed_mass_gap(Lp)
            assert gh > 0 and gp < 0
            print("曲率盒上正负权质量最小有向间隙：高度盒", gh.str(3), "相位盒", gp.str(3))
        s, bits = FAIL_HEIGHT[name]
        g = dyadic(gamma+arb(s), bits)
        assert test(at_height(g)) is False
        dh = abs(g-gamma)
        assert arb(2)**-kh < dh <= arb(2)**-(kh+HEIGHT_DROP[name])
        dx, dy = FAIL_PHASE[name]
        tx, ax = cayley_parameter(x0+arb(dx), 48)
        ty, ay = cayley_parameter(y0+arb(dy), 48)
        assert test(at_cayley(gamma, tx, ty)) is False
        dp = abs(ax-x0).max(abs(ay-y0))
        assert arb(2)**-kp < dp <= arb(2)**-(kp+PHASE_DROP[name])
        lo, hi = DIST[name]
        assert arb(lo[0]) < dh < arb(lo[1]) and arb(hi[0]) < dp < arb(hi[1])
        rows[name] = (kh, dh, kp, dp, g, tx, ty)
        print(name, "高度球 2^-%d 认证；失败高度距离" % kh, dh.str(8), "（2^%.2f）" % log2(dh),
              "| 相位球 2^-%d 认证；失败相位最大模距离" % kp, dp.str(8), "（2^%.2f）" % log2(dp))
    report["8359"]["rows"] = rows
    # 相邻零点
    for m in (8358, 8360):
        g = acb.zeta_zero(m).imag
        L = at_height(g)
        gap = phase_gap(gamma, g)
        vals = [test(L) for _, test in WITNESSES]
        report["8359"][m] = (g-gamma, gap, vals)
        # 门槛与相邻零点间距之比：三模（相位、高度）与曲率（相位）。
        assert gap/arb("0.05251") > 4.5 and gap/arb("0.00151") > 150
        assert abs(g-gamma)/arb("0.14") > 3.6
        (a, b), (c, d), expected = NEIGHBOUR[m]
        assert arb(a) < g-gamma < arb(b) and arb(c) < gap < arb(d) and vals == expected
        print("相邻第", m, "枚：高度差", (g-gamma).str(10), "相位差最大模", gap.str(10), "三模/曲率/强耦合", vals)


FLAG_NEIGHBOUR = {-1: (("-0.249296", "-0.249295"), ("0.0138799", "0.01388"), ("-0.12923", "-0.12922"), ("-0.082446", "-0.082445")),
                  1: (("0.0744657", "0.0744658"), ("-0.004147", "-0.004146"), ("0.04102", "0.041021"), ("0.024727", "0.024728"))}


def part_flagship(report):
    n = 2_306_888_100_085
    ctx.prec = 256
    gamma = acb.zeta_zero(n).imag
    flavor, rel, window, interval, K, off = flavor_tools()
    clocks = {d: acb(0, A[d]*gamma).exp() for d in A}
    x0, y0 = clocks[724].arg(), clocks[804].arg()
    delta0 = x0-y0
    c0 = flavor["singlet_at"](clocks)[1]
    assert c0.overlaps(K*(delta0-off).cos())
    u0, j0 = rel(c0)
    # 与原复算一致：+1.328335 ppm 与 +0.811085 ppm。
    assert arb("1.328334e-6") < u0 < arb("1.328336e-6") and arb("8.11084e-7") < j0 < arb("8.11085e-7")
    dA = A[724]-A[804]
    kh, kp = 20, 25
    assert interval(delta0, abs(dA)*arb(2)**-kh) is True
    assert interval(delta0, 2*arb(2)**-kp) is True
    # 同一方法在粗一位上不能认证。
    assert interval(delta0, abs(dA)*arb(2)**-(kh-1)) is not True
    assert interval(delta0, 2*arb(2)**-(kp-1)) is not True
    # 高度失败点：二进有理高度。
    g = dyadic(gamma+arb("1.3e-6"), 30)
    cg = flavor["singlet_at"]({d: acb(0, A[d]*g).exp() for d in A})[1]
    assert window(cg) is False
    dh = abs(g-gamma)
    # 相位失败点：两道钟各取有理 Cayley 点。
    tx, ax = cayley_parameter(x0-arb("3.5e-8"), 64)
    ty, ay = cayley_parameter(y0+arb("3.5e-8"), 64)
    cp_ = flavor["singlet_at"]({724: cayley(tx), 804: cayley(ty)})[1]
    assert window(cp_) is False
    dp = abs(ax-x0).max(abs(ay-y0))
    # 相位失败点的相对钟角偏移，以及高度失败点的相对钟角偏移。
    shift_p, shift_h = abs((ax-ay)-delta0), abs(dA)*dh
    assert shift_p < arb("7.01e-8") and arb("7.23e-8") < shift_h < arb("7.24e-8")
    assert arb(2)**-kh < dh < arb(2)**-(kh-1) and arb(2)**-kp < dp < arb(2)**-(kp-1)
    assert arb("1.29915e-6") < dh < arb("1.29916e-6") and arb("3.4999e-8") < dp < arb("3.5001e-8")
    ug, jg = rel(cg)
    up, jp = rel(cp_)
    assert arb("2.03426e-6") < ug < arb("2.03427e-6") and arb("1.24212e-6") < jg < arb("1.24213e-6")
    assert arb("2.01150e-6") < up < arb("2.01151e-6") and arb("1.22822e-6") < jp < arb("1.22823e-6")
    assert gamma.rad() < arb(2)**-215
    assert 2**39 <= gamma < 2**40 and arb("597041689632.5519") < gamma < arb("597041689632.552")
    print("旗舰味窗：γ", gamma.str(25), "原证书 γ 半径 2^", round(log2(gamma.rad()), 1))
    print("  高度球 2^-%d 认证，失败高度距离" % kh, dh.str(8), "（2^%.2f）" % log2(dh),
          "；相位球 2^-%d 认证，失败相位距离" % kp, dp.str(8), "（2^%.2f）" % log2(dp))
    print("  失败点读数 U、J 相对偏差", *[v.str(8) for v in rel(cg)], "|", *[v.str(8) for v in rel(cp_)])
    neighbours = {}
    for m in (n-1, n+1):
        gm = acb.zeta_zero(m).imag
        pm = {d: acb(0, A[d]*gm).exp() for d in A}
        cm = flavor["singlet_at"](pm)[1]
        dm = ((pm[724]/pm[804]).arg()-delta0)
        neighbours[m] = (gm-gamma, dm, rel(cm), window(cm))
        assert window(cm) is False
        assert abs(dm)/arb("7.01e-8") > arb("5.7e4")
        (a, b), (c, d), (e, f), (gg, hh) = FLAG_NEIGHBOUR[m-n]
        u, j = rel(cm)
        assert arb(a) < gm-gamma < arb(b) and arb(c) < dm < arb(d)
        assert arb(e) < u < arb(f) and arb(gg) < j < arb(hh)
        print("  相邻第", m, "枚：高度差", (gm-gamma).str(10), "相对钟角差", dm.str(10),
              "U、J 相对偏差", *[v.str(6) for v in rel(cm)])
    report["flag"] = dict(gamma=gamma, u=u0, j=j0, dh=dh, dp=dp, nb=neighbours)
    ctx.prec = 128


def main():
    t = time.time()
    report = {}
    part_8359(report)
    part_flagship(report)
    print("总用时 %.1f 秒" % (time.time()-t))


if __name__ == "__main__":
    main()
