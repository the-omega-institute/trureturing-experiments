"""二十七项字典、三枚零点 #11661、#11934、#22027 共用同一组系数：重—重上型块较小的精确共用完整命中。

区间证书，448 位。数据文件给出 72 个实坐标（w 的实部与虚部、三个真空 u_n 的实部与虚部），其中 33 个固定为有限小数，
其余 39 个构成 Krawczyk 盒。完整命中的 39 条实多项式方程（每零点 13 条）沿用
ch39-shared-zero-exact-hit-and-joint-heavy-bound.py 的 equations；Krawczyk 包含 K(X)⊂int X 认证盒内有唯一解，
该解是精确的三零点共用完整命中。随后在盒上包围九个质量、CKM 四量、重—重上型块（谱范数的 Rayleigh 下界与
Sylvester 上界、Frobenius 区间）、h² 系数、完整未投影上型顶角与重双重态×轻单态混合块的谱范数。
"""
import json
from fractions import Fraction
import sys
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

HERE = Path(__file__).parent
SH = run_path(str(HERE/"ch39-shared-zero-exact-hit-and-joint-heavy-bound.py"))
prev, single, MIXING = SH["prev"], SH["single"], SH["MIXING"]
adjoint, study = prev["adjoint"], prev["study"]
DATA = json.loads((HERE/"ch39-shared-zero-low-heavy-exact-hit.json").read_text())
ZEROS = (11661, 11934, 22027)

# 精确共用命中的认证区间：重—重块谱范数（Rayleigh 下界, Sylvester 上界）、Frobenius、h² 系数 GeV²，
# 以及完整未投影上型顶角与重双重态×轻单态混合块谱范数的 Rayleigh 下界。
CLAIMS = {
    "wnorm": ("7102.87", "7102.88"),
    "11661": dict(spec=("9.7624", "9.7625"), frob=("19.9970", "19.9971"), h2=("3.163e33", "3.165e33"), full="83.0", mixed="82.5"),
    "11934": dict(spec=("9.7235", "9.7236"), frob=("19.8669", "19.8670"), h2=("3.147e33", "3.149e33"), full="83.4", mixed="82.9"),
    "22027": dict(spec=("9.8161", "9.8162"), frob=("20.2056", "20.2057"), h2=("3.185e33", "3.187e33"), full="84.0", mixed="83.6"),
    "tau_upper": "9.8162",
}
OLD_SHARED_SPEC_LO = "65983.2"     # 已有精确共用命中三个块谱范数下界中最小者（引自其脚本的 HIT 表）
OLD_SHARED_H2_LO = "5.748e40"      # 已有精确共用命中 h² 系数下界中最小者
PERZERO_H2_FLOOR = "1.130e31"     # 逐零点证书 h² 系数下界的最小者（引自 ch39-per-leg-double-hop-transfer.json）
SINGLE_SPEC_HI = ("4.4319", "4.4289", "4.1574")   # 单零点精确命中的重块谱范数上界（引自单零点脚本的 CLAIMS）


def krawczyk(A0s):
    x = [arb(s) for s in DATA["x"]]
    active = DATA["active"]
    assert len(x) == 72 and len(active) == 39 and len(set(active)) == 39 and all(0 <= a < 72 for a in active)
    rad = arb(DATA["radius"])
    vals, grads = SH["equations"](A0s, x, active)
    F = arb_mat([[v] for v in vals])
    Jm = arb_mat([[arb(g.mid()) for g in row] for row in grads])
    one = arb_mat([[int(i==j) for j in range(39)] for i in range(39)])
    Y = Jm.solve(one, algorithm="approx")
    Y = arb_mat([[arb(Y[i,j].mid()) for j in range(39)] for i in range(39)])
    box, half = list(x), []
    for a in active:
        r = rad*arb(max(arb(1), abs(x[a])).upper())
        half.append(arb(0, r.upper()))
        box[a] = x[a] + half[-1]
    _, gb = SH["equations"](A0s, box, active)
    K = -(Y*F) + (one - Y*arb_mat(gb))*arb_mat([[box[a]-x[a]] for a in active])
    assert all(abs(K[i,0]).upper() < half[i].rad() for i in range(39))
    ratio = max((abs(K[i,0])/half[i].rad()).upper() for i in range(39))
    assert ratio < arb("3e-4")
    print("Krawczyk：39 维盒（半宽", DATA["radius"], "·max(1,|x_a|)）内 K(X)⊂int X，分量比最大", arb(ratio).str(3))
    return box


def spectral_interval(Y, top, hi):
    # 下端：Rayleigh 商（给定向量）；上端：hi²·1 − Y†Y ≻ 0（Sylvester）。
    single["spectral_below"](Y, hi)
    return single["rayleigh"](Y, top)


def measure(models, box):
    T = single["target_masses"]()
    target = [y for g in ("up", "down", "leptons") for y in T[g]]
    goals = [arb(a).sin() for a in (".22704", ".04275", ".00377")] + [arb("1.139")]
    w = acb_mat([[acb(box[t], box[27+t])] for t in range(27)])
    out = {"wnorm": sum((abs(w[t,0])**2 for t in range(27)), arb(0)).sqrt()}
    for a, n in enumerate(ZEROS):
        o = 54+6*a
        u = acb_mat([[acb(box[o+k], box[o+3+k])] for k in range(3)])
        z = u/(adjoint(u)*u)[0,0].real.sqrt()
        blocks, heavy, cols = models[n]
        matrices = study["mass_matrices"](blocks, w, z)
        masses = []
        for y in matrices:
            sv = sorted((r.real.sqrt() for r in (y*adjoint(y)).eig(algorithm="rump")), key=float)
            assert all(sv[i].upper() < sv[i+1].lower() for i in range(2))
            masses += sv
        angles = MIXING["observables"](*matrices[:2])[4:]
        # 盒内的精确解使每个包围都含其精确值；方程给出质量与 CKM 模方精确等于目标，相位方程给出 δ≡1.139 (mod π)；
        # δ 的包围含 1.139 而不含 1.139−π，于是 δ 精确等于 1.139。
        assert all((y-ref).contains(0) for y, ref in zip(masses, target))
        assert all((v-ref).contains(0) for v, ref in zip(angles[:4], goals))
        assert not (angles[3] - arb("1.139") + arb.pi()).contains(0)
        width = max(abs(y/ref-1).upper() for y, ref in zip(masses, target))
        mix = max(abs(v-ref).upper() for v, ref in zip(angles[:4], goals))
        assert width < arb("1e-21") and mix < arb("1e-18")
        assert abs(angles[4] - arb("3.20630e-5")) < arb("5e-11")
        flat = sum((z[k,0]*(cols[k]*w) for k in range(3)), acb_mat(361, 1))
        Yh = acb_mat([[flat[19*s+t,0] for t in range(19)] for s in range(19)])
        full, hh, mixed = single["vertex_blocks"](n, w, z)
        assert all((hh[i,j]-Yh[i,j]).contains(0) for i in range(19) for j in range(19))
        frob2 = sum((abs(Yh[i,j])**2 for i in range(19) for j in range(19)), arb(0))
        (q_lo, q_hi), (u_lo, u_hi) = single["heavy_mass_balls"](n)
        quad = single["quad"]
        out[n] = dict(Yh=Yh, full=full, mixed=mixed, frob=frob2.sqrt(),
                      h2=(quad(q_lo**2, u_lo**2)*frob2, quad(q_hi**2, u_hi**2)*frob2),
                      width=width, mix=mix, J=angles[4], z=[abs(z[k,0]) for k in range(3)])
    return out


def top_vector(Y):
    # 由中心矩阵的数值奇异向量给出 Rayleigh 商用的向量（只影响下界的紧度，不影响其正确性）。
    import numpy as np
    m = np.array([[complex(float(Y[i,j].real.mid()), float(Y[i,j].imag.mid())) for j in range(19)] for i in range(19)])
    v = np.linalg.svd(m)[2][0].conj()
    return [[repr(float(c.real)), repr(float(c.imag))] for c in v]


def verify(discover=False):
    ctx.prec = 448
    models = {n: prev["model"](n) for n in ZEROS}
    A0s = {n: SH["coefficient_matrix"](models[n][0]) for n in ZEROS}
    box = krawczyk(A0s)
    meas = measure(models, box)
    four_pi = 4*arb.pi()
    if discover:
        print("wnorm", meas["wnorm"].str(12))
        for n in ZEROS:
            m = meas[n]
            for key in ("Yh", "full", "mixed"):
                print(n, key, single["rayleigh"](m[key], top_vector(m[key])).str(12))
            print(n, "frob", m["frob"].str(12), "h2", m["h2"][0].str(6), m["h2"][1].str(6), "width", m["width"].str(3), m["mix"].str(3), "J", m["J"].str(10))
        return
    lo_w, hi_w = map(arb, CLAIMS["wnorm"])
    assert lo_w < meas["wnorm"] < hi_w
    for n in ZEROS:
        m, c = meas[n], CLAIMS[str(n)]
        lo, hi = map(arb, c["spec"])
        assert spectral_interval(m["Yh"], top_vector(m["Yh"]), c["spec"][1]) > lo
        assert arb(c["frob"][0]) < m["frob"] < arb(c["frob"][1])
        assert arb(c["h2"][0]) < m["h2"][0] and m["h2"][1] < arb(c["h2"][1])
        assert single["rayleigh"](m["full"], top_vector(m["full"])) > arb(c["full"])
        assert single["rayleigh"](m["mixed"], top_vector(m["mixed"])) > arb(c["mixed"])
        assert arb(c["full"]) > four_pi and arb(c["mixed"]) > four_pi
        print(n, "九质量与 CKM 四量包围均含目标（质量包围与目标的相对偏差 <", m["width"].str(3), "，CKM 包围与目标的偏差 <", m["mix"].str(3), "），J", m["J"].str(8))
        print("  重—重块谱范数 ∈ (", c["spec"][0], ",", c["spec"][1], ")，Frobenius ∈ (", c["frob"][0], ",", c["frob"][1],
              ")，h² 系数 ∈ (", c["h2"][0], ",", c["h2"][1], ") GeV²")
        print("  完整未投影上型顶角谱范数 >", c["full"], "，重双重态×轻单态块谱范数 >", c["mixed"])
    top = arb(max(arb(CLAIMS[str(n)]["spec"][1]).upper() for n in ZEROS))      # 十进制端点的最大者
    bottom = arb(min(arb(CLAIMS[str(n)]["spec"][0]).lower() for n in ZEROS))
    # τ0 的上界取十进制上端的最大者（以有理数比较，不比较球）。
    assert max(Fraction(CLAIMS[str(n)]["spec"][1]) for n in ZEROS) == Fraction(CLAIMS["tau_upper"])
    return meas, top, bottom, four_pi


def comparisons(meas, top, bottom, four_pi):
    # 三个块谱范数的上端都小于 4π；下端都大于 4.44 与 √(4π)。所以 τ0<4π 已判定，τ0 与 4.44、√(4π) 的大小关系未判定。
    assert top < four_pi
    assert bottom > arb("4.44") and bottom > four_pi.sqrt()
    # 与已有结果的核对：引用的常数逐项与其脚本中的表一致。
    assert SH["HIT"]["22027"]["spec"][0] == OLD_SHARED_SPEC_LO and SH["HIT"]["22027"]["quad"][0] == OLD_SHARED_H2_LO
    assert min(Fraction(SH["HIT"][str(n)]["quad"][0]) for n in ZEROS) == Fraction(OLD_SHARED_H2_LO)
    # 逐零点证书给出的任何命中的 h² 系数下界，三者最小者为 1.130e31（#11934），引自其数据文件，本脚本不复算。
    floors = prev["DATA"]["all_bounds"]
    assert min(Fraction(floors[str(n)][2]) for n in ZEROS) == Fraction(PERZERO_H2_FLOOR) == Fraction(floors["11934"][2])
    assert all(Fraction(SH["HIT"][str(n)]["spec"][0]) >= Fraction(OLD_SHARED_SPEC_LO) for n in ZEROS)
    assert tuple(single["CLAIMS"][str(n)]["spec"][1] for n in ZEROS) == SINGLE_SPEC_HI
    lower = SH["JOINT"]["11661,11934,22027"][2]          # 联合证书下界 0.6233（引自已有脚本，本脚本不复算）
    assert lower == ".6233"
    # 两零点集合：联合证书下界（引自已有脚本）与本命中给出的上界（块谱范数十进制上端的最大者，以有理数比较）。
    pairs = ("11661,11934", "11661,22027", "11934,22027")
    assert tuple(SH["JOINT"][k][2] for k in pairs) == (".5742", ".5944", ".5342")
    for k, ub in zip(pairs, ("9.7625", "9.8162", "9.8162")):
        assert max(Fraction(CLAIMS[n]["spec"][1]) for n in k.split(",")) == Fraction(ub)
    # 已有精确共用命中的真空近于 (0,0,1)。
    ox = [arb(v) for v in SH["DATA"]["hit"]["x"]]
    for a in range(3):
        o = 54+6*a
        mod = [(ox[o+k]**2 + ox[o+3+k]**2).sqrt() for k in range(3)]
        r = (mod[0]**2 + mod[1]**2 + mod[2]**2).sqrt()
        assert mod[0]/r < arb("1e-7") and mod[1]/r < arb("1e-6") and mod[2]/r > 1 - arb("1e-12")
    ratio = arb(OLD_SHARED_SPEC_LO)/top
    assert ratio > 6700
    assert SH["WNORM"][0] == "6795043" and arb(SH["WNORM"][0])/arb(CLAIMS["wnorm"][1]) > 950
    # 真空分量模：三枚零点都接近纯特征标真空 (0,1,0)。
    for n in ZEROS:
        a0, a1, a2 = meas[n]["z"]
        assert a0 < arb(".03") and a1 > arb(".9997") and a2 < arb(".007")
    h2_hi = max(arb(CLAIMS[str(n)]["h2"][1]).upper() for n in ZEROS)
    h2_lo = min(arb(CLAIMS[str(n)]["h2"][0]).lower() for n in ZEROS)
    assert arb(OLD_SHARED_H2_LO)/arb(h2_hi) > arb("1.8e7")
    assert arb(h2_lo)/arb("125.25")**2 > arb("2.0e29") and arb(h2_lo) > arb(PERZERO_H2_FLOOR)
    single_top = arb(max(arb(x).upper() for x in SINGLE_SPEC_HI))
    assert top/single_top < arb("2.22")
    print("τ0 = 精确共用命中上 max_n‖Y_hh^(n)‖₂ 的下确界 ∈ (", lower, ",", CLAIMS["tau_upper"], ")；4π =", four_pi.str(6))
    print("  该命中三个块均 > 4.44 > √(4π)；与已有精确共用命中相比上界缩小 >", 6700, "倍；h² 系数 <", arb(h2_hi).str(4), "GeV²")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--discover":
        verify(discover=True)
    else:
        comparisons(*verify())
