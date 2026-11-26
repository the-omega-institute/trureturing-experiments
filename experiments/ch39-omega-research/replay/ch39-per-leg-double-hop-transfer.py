"""每条腿至多两次邻边输运（新作用律）时，三枚候选共用系数的质量比障碍、完整命中见证与重块下界。"""
import json
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


HERE = Path(__file__).parent
study = run_path(str(HERE/"ch39-three-leg-yukawa.py"))
edge, shell = study["edge"], study["shell"]
adjoint, positive = study["adjoint"], study["positive"]
DATA = json.loads((HERE/"ch39-per-leg-double-hop-transfer.json").read_text())

ZEROS = (11661, 11934, 22027)
MONOS = tuple(product(range(3), repeat=3))          # 二十七项字典：每腿 T^0,T^1,T^2
EIGHT = tuple(product(range(2), repeat=3))          # 原八项
EXTRA = tuple(m for m in MONOS if m not in EIGHT)   # 十九个含二次输运的单项
RD, RE = arb("61.36519"), arb("84.98147")


def legs(number):
    # 与 ch39-three-leg-yukawa.py 相同的低模、相位约定与输运 T=-offdiag(L_q)/8。
    height = acb.zeta_zero(number).imag
    neutral = shell["laplacian"](0, height)
    center = sum(neutral[i,i].real for i in range(19)) / 18
    fields, ops = {}, {}
    one = acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    for q in (1, -4, 2, -3, 6, 3):
        lap = shell["laplacian"](q, height)
        T = acb_mat([[-lap[i,j]/8 if i != j else 0 for j in range(19)] for i in range(19)])
        modes = edge["mode_columns"](edge["filters"]["spectra"](q, height, center))
        for j in range(3):
            anchor = next(modes[i,j] for i in range(19) if abs(modes[i,j]) > arb(".01"))
            phase = anchor.conjugate()/abs(anchor)
            for i in range(19):
                modes[i,j] *= phase
        fields[q] = [modes, T*modes, T*(T*modes)]
        ops[q] = [one, T, T*T]
    return fields, ops


def sector(left, right, higgs, conjugate):
    rows = [[] for _ in range(3)]
    for k in range(3):
        for i in range(3):
            j = (-i + (k if conjugate else -k)) % 3
            rows[k].append([sum(left[a][x,i]*right[b][x,j]*higgs[c][x,k] for x in range(19))
                            for a,b,c in MONOS])
    return [acb_mat(r) for r in rows]


def heavy_grams(fields, ops):
    # 每个 Higgs 特征标 k 上，重双重态×重单态上型块对二十七个系数的 Gram 矩阵。
    h, vq, vu = fields[3], fields[1][0], fields[-4][0]
    one = ops[1][0]
    left = one - vq.conjugate()*vq.transpose()
    right = one - vu*adjoint(vu)
    grams, columns = [], []
    for k in range(3):
        cols = []
        for a,b,c in MONOS:
            col = ops[1][a].transpose()*acb_mat([[h[c][x,k]*ops[-4][b][x,j] for j in range(19)]
                                                  for x in range(19)])
            cols.append(list((left*col*right).entries()))
        C = acb_mat(cols).transpose()
        columns.append(C)
        grams.append(adjoint(C)*C)
    return grams, columns


def heavy_masses(number):
    # 最轻重双重态与重单态质量，单位与 ch39-three-leg-yukawa.py 相同。
    height = acb.zeta_zero(number).imag
    neutral = shell["laplacian"](0, height)
    center = sum(neutral[i,i].real for i in range(19)) / 18
    out = []
    for q in (1, -4):
        levels = [abs(r.real*(r.real-center)) for r in shell["laplacian"](q, height).eig(algorithm="rump")]
        assert all(v > 9 or v < 9 for v in levels)
        heavy = [v for v in levels if v > 9]
        assert len(heavy) == 16
        out.append(min(v.lower() for v in heavy)*arb("2.435e18")/724**2*(32 if q == -4 else 1))
    return out


def lesser(a, b):
    # 两球之最小值的包络：下端取两下端之小，上端取两上端之小，不依赖球比较。
    lo, hi = min(a.lower(), b.lower()), min(a.upper(), b.upper())
    return arb(lo).union(arb(hi))


def sub(matrix, idx, cols_only=False):
    rows = range(matrix.nrows()) if cols_only else idx
    return acb_mat([[matrix[i,j] for j in idx] for i in rows])


def column(values):
    return acb_mat([[acb(arb(x), arb(y))] for x,y in values])


def model(number):
    fields, ops = legs(number)
    hbar = [b.conjugate() for b in fields[3]]
    blocks = list(zip(sector(fields[1], fields[-4], fields[3], False),
                      sector(fields[1], fields[2], hbar, True),
                      sector(fields[-3], fields[6], hbar, True)))
    heavy, heavy_cols = heavy_grams(fields, ops)
    return blocks, heavy, heavy_cols


def nine_term_witnesses(models):
    # 原八项加任一含二次输运的单项，纯特征标 2 真空，六个范数比逐一认证。
    reach_one = []
    for key, entry in DATA["nine"].items():
        e = tuple(json.loads(key))
        idx = [MONOS.index(m) for m in EIGHT+(e,)]
        w = column(entry["w"])
        worst = None
        for number in ZEROS:
            u, d, l = [(adjoint(sub(c, idx, True)*w)*(sub(c, idx, True)*w))[0,0].real
                       for c in models[number][0][2]]
            for ratio in ((u/d).sqrt()/RD, (u/l).sqrt()/RE):
                worst = ratio if worst is None else lesser(worst, ratio)
        assert worst > arb(".987")
        if key in DATA["nine_values"]:
            lo, hi = map(arb, DATA["nine_values"][key])
            assert lo < worst < hi
        if worst > 1:
            reach_one.append(e)
        print("九项", e, "六个比与参考下限之比的最小值", worst.str(8))
    assert len(reach_one) == 15
    return reach_one


def leakage(models):
    # 不依赖正交性的严格修正：‖Y‖²≥(1−2ε/λ)Σp_k w†H_k w，ε 为交叉 Gram 的 Frobenius 上界，λ 为 H_k 的下界。
    out = {}
    for number in ZEROS:
        _, heavy, cols = models[number]
        eps = arb(0)
        for k in range(3):
            for l in range(k+1, 3):
                cross = adjoint(cols[k])*cols[l]
                square = sum((abs(cross[i,j])**2 for i in range(27) for j in range(27)), arb(0))
                frob = arb(square.upper()).sqrt()          # 平方和的上端点，再开方
                assert frob.is_finite()
                eps = arb(max(eps.upper(), frob.upper()))   # 端点比较取上界的最大者
        lam = arb(DATA["heavy_floor"][str(number)])
        one = acb_mat([[int(i==j) for j in range(27)] for i in range(27)])
        for k in range(3):
            positive(heavy[k] - lam*one)
        delta = 2*eps/lam
        assert eps < arb("9.95e-51") and delta < arb("5.81e-42")
        print(number, "交叉 Gram Frobenius 上界", eps.str(3), "H_k 下界", lam.str(3), "修正因子 1−", delta.str(3))
        out[number] = 1 - delta
    return out


def caps(models, minima, factor):
    # cap·H_k−U_k+μ(D_k−U_k)≻0 对 k=0,1,2 共用一个 μ≥0。
    refs = run_path(str(HERE/"ch39-shell-flavor-matching-windows.py"))
    assert refs["up"][-1]*arb(".999") > arb(".86")   # 命中者的上型 Frobenius 范数 ≥ m_t > 0.86
    results = {}
    for key, row in DATA["caps"].items():
        idx = list(range(27)) if key == "all" else [MONOS.index(m) for m in EIGHT+(tuple(json.loads(key)),)]
        worst = None
        for number in ZEROS:
            cap, mu = map(arb, row[str(number)])
            assert cap > 0 and mu >= 0
            blocks, heavy, _ = models[number]
            for k in range(3):
                u, d, l = [sub(adjoint(c)*c, idx) for c in blocks[k]]
                pair = (RD**2*d + RE**2*l)/2
                positive(cap*sub(heavy[k], idx) - u + mu*(pair - u))
            frob2 = factor[number]*arb(".86")**2/cap
            spectral = (frob2/16).sqrt()
            worst = spectral if worst is None else lesser(worst, spectral)
            if key == "all":
                q2, u2 = [m**2 for m in minima[number]]
                quadratic = (q2+u2-q2*u2/(q2+u2))*frob2
                f_lo, s_lo, h_lo = map(arb, DATA["all_bounds"][str(number)])
                assert frob2 > f_lo and spectral > s_lo and quadratic > h_lo
                assert quadratic/arb("125.25")**2 > arb("7e26")
                print("二十七项", number, "重块 Frobenius² 下界", frob2.str(6),
                      "谱范数下界", spectral.str(6), "四次迹 h² 系数下界 GeV²", quadratic.str(6))
        results[key] = worst
        if key != "all":
            print("九项", tuple(json.loads(key)), "三枚零点中最弱的强制谱范数下界", worst.str(6))
    return results


def rotation_consistency():
    # 相容性核对（非证明）：坐标轮换置换 R 与每个 T_q 对易的包络含零。重块下界不依赖它，见 leakage。
    sites = shell["geo"]["sites"]
    rotate = [sites.index((b,c,a)) for a,b,c in sites]
    for number in ZEROS:
        height = acb.zeta_zero(number).imag
        for q in (1, -4, 3):
            lap = shell["laplacian"](q, height)
            assert all((lap[rotate[i],rotate[j]]-lap[i,j]).contains(0) for i in range(19) for j in range(19))


def full_witness(models, minima):
    refs = run_path(str(HERE/"ch39-shell-flavor-matching-windows.py"))
    assert refs["up"][-1]*arb(".999") > arb(".86")
    mixing = run_path(str(HERE/"ch39-finite-orthogonal-ckm.py"))
    target = [y for g in ("up","down","leptons") for y in refs[g]]
    goals = [arb(x).sin() for x in (".22704",".04275",".00377")] + [arb("1.139")]
    w = column(DATA["full"]["w"])
    for number in ZEROS:
        z = column(DATA["full"]["higgs"][str(number)])
        z = z/(adjoint(z)*z)[0,0].real.sqrt()
        blocks, heavy, _ = models[number]
        matrices = study["mass_matrices"](blocks, w, z)
        masses = [s for y in matrices for s in sorted(
            (r.real.sqrt() for r in (y*adjoint(y)).eig(algorithm="rump")), key=float)]
        error = max(abs(y/ref-1).upper() for y,ref in zip(masses, target))
        angles = mixing["observables"](*matrices[:2])[4:]
        mix = max(abs(v-ref).upper() for v,ref in zip(angles[:4], goals))
        m_tol, x_tol = map(arb, DATA["full"]["tol"][str(number)])
        assert error < m_tol and mix < x_tol
        assert abs(angles[4] - arb("3.20630e-5")) < arb("5e-11")
        # 见证本身的重—重上型块 Y=Σ_k z_k Σ_t w_t Y_{k,t}，直接组装。
        cols = models[number][2]
        flat = sum((z[k,0]*(cols[k]*w) for k in range(3)), acb_mat(361,1))
        Y = acb_mat([[flat[19*s+t,0] for t in range(19)] for s in range(19)])
        frob2 = sum((abs(Y[i,j])**2 for i in range(19) for j in range(19)), arb(0))
        gram = adjoint(Y)*Y
        v = column(DATA["full"]["top"][str(number)])
        rayleigh = ((adjoint(v)*gram*v)[0,0].real/(adjoint(v)*v)[0,0].real).sqrt()
        q2, u2 = [m**2 for m in minima[number]]
        quadratic = (q2+u2-q2*u2/(q2+u2))*frob2
        f_lo, s_lo, h_lo = map(arb, DATA["full"]["heavy"][str(number)])
        assert frob2.sqrt() > f_lo and rayleigh > s_lo and quadratic > h_lo
        assert rayleigh > frob2.sqrt()/4
        print("二十七项共同系数", number, "九质量最大相对误差", error.str(4),
              "CKM 四量最大绝对误差", mix.str(4), "J", angles[4].str(6))
        print("  见证重块 Frobenius 范数", frob2.sqrt().str(6), "谱范数下界（Rayleigh）", rayleigh.str(6),
              "按原式的 h² 系数下界 GeV²", quadratic.str(6))


def verify():
    ctx.prec = 448
    models = {n: model(n) for n in ZEROS}
    rotation_consistency()
    reach = nine_term_witnesses(models)
    print("达到全部参考下限的九项扩充数", len(reach))
    minima = {n: heavy_masses(n) for n in ZEROS}
    for n in ZEROS:
        assert arb("5.08e13") < minima[n][0] < arb("5.10e13") and arb("2.79e15") < minima[n][1] < arb("2.83e15")
        print(n, "最轻重双重态与重单态质量 GeV", *(m.str(5) for m in minima[n]))
    factor = leakage(models)
    bounds = caps(models, minima, factor)
    nine = {k: v for k,v in bounds.items() if k != "all"}
    low = None
    for v in nine.values():
        low = v if low is None else lesser(low, v)
    assert low > arb("2.7")
    above = [k for k,v in nine.items() if v > arb(4*arb.pi()).sqrt()]
    assert len(above) == 18 and "[1, 2, 2]" not in above
    special = {"[1, 2, 2]": "2.719", "[2, 2, 2]": "5.79", "[2, 2, 0]": "7.08", "[0, 2, 1]": "13.04",
               "[2, 0, 0]": "14.12", "[2, 0, 1]": "14.21", "[2, 1, 1]": "10.23"}
    for k, v in special.items():
        assert nine[k] > arb(v)
    others = [v for k,v in nine.items() if k not in special]
    assert all(arb("9.25") < v < arb("16.01") for v in others)
    print("十九个九项扩充在每枚零点上的重块谱范数下界的最小值", low.str(6),
          "超过 √(4π) 的扩充数", len(above))
    print("二十七项字典在三枚零点上的谱范数下界", bounds["all"].str(6))
    full_witness(models, minima)


if __name__ == "__main__":
    verify()
