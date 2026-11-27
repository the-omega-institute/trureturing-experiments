"""五荷两步三维轻空间在整张联合相位环面上的两两泄漏：三重轮换分解与 Arb 盒证书。

记号见同名说明。荷 q 在联合相位 (X, Y) 上读 (e^{iqX}, e^{iqY})；对辐边亏额
delta 取遍 [0, delta_max]，delta_max = (23/25) 2^-20 / (14^2 + 30)。
"""
import sys
import time
from fractions import Fraction as Fr
from itertools import combinations
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, ctx

ROOT = Path(__file__).resolve().parent
SHELL = run_path(str(ROOT / "ch39-three-ray-chiral-boundary.py"))
GEO, SOFT = SHELL["geo"], SHELL["soft"]
ctx.prec = 96

SITES, EDGES, POLE, TUNED = GEO["sites"], GEO["edges"], SOFT["pole"], SOFT["star"]["orbit"]
CHARGES = (1, 2, -3, -4, 6)
ABS = (1, 2, 3, 4, 6)
DMAX = arb(23)/25*arb(2)**-20/(14**2+30)
CBALL = arb(15)-DMAX/6+arb(0, 1)*DMAX/6          # c = 15 - delta/3 的全体取值
LIGHT = {s: {1: [2], 2: [2], 0: ([3] if s % 2 else [0, 3])} for s in ABS}


# ---------------- 三重轮换与扇区基 ----------------
ROT = [SITES.index((y, z, x)) for x, y, z in SITES]
ORBITS, _seen = [], set()
for _i in range(19):
    if _i not in _seen:
        _o = [_i]
        while ROT[_o[-1]] != _i:
            _o.append(ROT[_o[-1]])
        _seen.update(_o)
        ORBITS.append(_o)
assert sorted(len(o) for o in ORBITS) == [1]+[3]*6
CENTER = [n for n, o in enumerate(ORBITS) if 6 in o][0]
assert ORBITS[CENTER] == [6, 7, 8] and [POLE] in ORBITS and tuple(ORBITS[0]) == tuple(TUNED)
OMEGA = acb(0, 2*arb.pi()/3).exp()


def sector_basis(k):
    cols = []
    for o in ORBITS:
        if len(o) == 1:
            if k == 0:
                cols.append({o[0]: acb(1)})
        else:
            cols.append({o[m]: OMEGA**(-k*m)/arb(3).sqrt() for m in range(3)})
    B = acb_mat(19, len(cols))
    for c, col in enumerate(cols):
        for i, v in col.items():
            B[i, c] = v
    return B, [n for n, o in enumerate(ORBITS) if k == 0 or len(o) == 3]


def dagger(M):
    return M.conjugate().transpose()


def parts(s):
    """L_s(U,V) = F + U C1 + conj(U) C1^H + V C2 + conj(V) C2^H - delta K。"""
    F, C1, C2, K = (acb_mat(19, 19) for _ in range(4))
    for i, j, d, z in EDGES:
        spoke = POLE in (i, j)
        for M in (F,):
            M[i, i] += 1
            M[j, j] += 1
        if spoke and (i in TUNED or j in TUNED):
            K[i, i] += 1
            K[j, j] += 1
            K[i, j] -= 1
            K[j, i] -= 1
        phase = acb(1) if spoke else z**s
        if d in (724, 804):
            assert i < 6 and j in (6, 7, 8) and not spoke
            (C1 if d == 724 else C2)[i, j] -= phase
        else:
            F[i, j] -= phase
            F[j, i] -= phase.conjugate()
    return F, C1, C2, K


BASES = [sector_basis(k) for k in range(3)]
ROTM = acb_mat([[int(ROT[i] == j) for i in range(19)] for j in range(19)])


def block_data(s):
    F, C1, C2, K = parts(s)
    for M in (F, C1, C2, K):
        assert (M*ROTM-ROTM*M).contains(acb_mat(19, 19))          # 与轮换对易
    out = []
    for k, (B, orbs) in enumerate(BASES):
        n = len(orbs)
        c = orbs.index(CENTER)
        Fk = dagger(B)*F*B
        g = []
        for C in (C1, C2):
            G = dagger(B)*C*B
            for a in range(n):
                for b in range(n):
                    if b != c:
                        assert G[a, b].contains(0) and abs(G[a, b]).upper() < arb(2)**-60
            assert G[c, c].contains(0)
            g.append([G[a, c] for a in range(n)])
        Kk = dagger(B)*K*B
        assert (Fk-dagger(Fk)).contains(acb_mat(n, n))
        norms = [sum((abs(x)**2 for x in v), arb(0)).sqrt() for v in g]
        out.append(dict(F=Fk, g=g, gn=norms, K=Kk, n=n, c=c,
                        Ff=np.array([[complex(Fk[a, b].mid()) for b in range(n)] for a in range(n)]),
                        gf=[np.array([complex(x.mid()) for x in v]) for v in g]))
    return out


BLOCKS = {s: block_data(s) for s in ABS}
# 软化星图 K（三条辐边 14-{0,3,4} 的拉普拉斯）：限制到 {14,0,3,4} 的首一特征多项式在五个整数点上
# 等于 x (x-1)^2 (x-4)，故谱为 {0,1,1,4}，||K|| = 4。
def _det(m):
    m = [row[:] for row in m]
    n, d = len(m), Fr(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if m[r][i] != 0), None)
        if piv is None:
            return Fr(0)
        if piv != i:
            m[i], m[piv], d = m[piv], m[i], -d
        d *= m[i][i]
        for r in range(i+1, n):
            f = m[r][i]/m[i][i]
            m[r] = [a-f*b for a, b in zip(m[r], m[i])]
    return d


_star = [POLE, *TUNED]
_K = parts(1)[3]
_K4 = [[Fr(round(float(_K[i, j].real.mid()))) for j in _star] for i in _star]
assert all(_K[i, j].contains(int(_K4[a][b])) for a, i in enumerate(_star) for b, j in enumerate(_star))
assert all(_K[i, j] == 0 for i in range(19) for j in range(19) if i not in _star or j not in _star)
assert all(_det([[Fr(x)*(a == b)-_K4[a][b] for b in range(4)] for a in range(4)]) == x*(x-1)**2*(x-4)
           for x in range(5))
assert all(abs(BLOCKS[1][k]["K"][a, b]-BLOCKS[1][k]["K"][b, a].conjugate()).upper() < 1e-20
           for k in range(3) for a in range(BLOCKS[1][k]["n"]) for b in range(BLOCKS[1][k]["n"]))


def center_matrix(blk, U, V):
    n, c = blk["n"], blk["c"]
    M = acb_mat(blk["F"])
    for vec, Z in zip(blk["g"], (U, V)):
        for a in range(n):
            M[a, c] += Z*vec[a]
            M[c, a] += (Z*vec[a]).conjugate()
    return M


def hval(lo, hi):
    lam = (lo+hi)/2+arb(0, 1)*(hi-lo)/2
    return abs(lam*(lam-CBALL))


def turn(t):
    return acb(0, 2*arb.pi()*arb(t.numerator)/t.denominator).exp()


def radians(h):
    return 2*arb.pi()*arb(h.numerator)/h.denominator


def eigen_sector(blk, U, V):
    """中心矩阵的近似本征对与残差。中心矩阵取 delta = 0。"""
    n, c = blk["n"], blk["c"]
    Uf, Vf = complex(U.mid()), complex(V.mid())
    Mf = blk["Ff"].copy()
    for vec, Z in zip(blk["gf"], (Uf, Vf)):
        Mf[:, c] += Z*vec
        Mf[c, :] += np.conj(Z*vec)
    mu, vecs = np.linalg.eigh((Mf+Mf.conj().T)/2)
    M = center_matrix(blk, U, V)
    Vm = acb_mat([[complex(vecs[a, b]) for b in range(n)] for a in range(n)])
    RES = M*Vm
    rho, nv, cols = [], [], []
    for j in range(n):
        col = [Vm[a, j] for a in range(n)]
        norm = sum((abs(x)**2 for x in col), arb(0)).sqrt()
        r = sum((abs(RES[a, j]-arb(mu[j])*col[a])**2 for a in range(n)), arb(0)).sqrt()/norm
        rho.append(r.upper())
        nv.append(norm)
        cols.append(col)
    # 埃尔米特残差定理：每个区间 [mu_j - rho_j, mu_j + rho_j] 含一个本征值；区间两两不交，
    # 故恰各含一个，且次序就是本征值的编号。
    for j in range(n-1):
        assert arb(mu[j])+rho[j] < arb(mu[j+1])-rho[j+1]
    return dict(mu=[arb(x) for x in mu], rho=rho, nv=nv, cols=cols)


def analyze(s, U, V, rad):
    """荷 s 在自身相位 (U, V) 的盒 {|dU|, |dV| <= rad}、全部 delta 上的读数。U, V 可为球，残差对球中每点成立。

    返回：
      sel   盒上经认证的两步轻线编号（每扇区一个；偶荷扇区 0 可能两个候选），失败为 None；
      lines 指定候选编号的近似向量及其与盒内任一点真本征线的角度上界；
      psel, plines 同上，但只对盒心（仍对全部 delta）。
    """
    secs = [eigen_sector(blk, U, V) for blk in BLOCKS[s]]
    out = {}
    for tag, r in (("box", rad), ("pt", arb(0))):
        H, lines = {}, {}
        for k, (blk, sec) in enumerate(zip(BLOCKS[s], secs)):
            n, c = blk["n"], blk["c"]
            eps = r*(blk["gn"][0]+blk["gn"][1])+4*DMAX          # ||M_p - M_c|| 的上界（Weyl）
            for j in range(n):
                H[k, j] = hval(sec["mu"][j]-sec["rho"][j]-eps, sec["mu"][j]+sec["rho"][j]+eps)
            for i in DESIGNATED[s][k]:
                col, nv = sec["cols"][i], sec["nv"][i]
                # ||E v|| <= |dU| (|v_c| ||g1|| + sum|g1_a v_a|) + |dV| (...) + 4 delta ||v||
                ev = abs(col[c])*(blk["gn"][0]+blk["gn"][1])
                for g in blk["g"]:
                    ev += sum((abs(g[a])*abs(col[a]) for a in range(n)), arb(0))
                res = arb(sec["rho"][i])+(r*ev+4*DMAX*nv)/nv
                sep = None
                for j in range(n):
                    if j != i:
                        d = abs(sec["mu"][i]-sec["mu"][j])-sec["rho"][j]-eps
                        sep = d if sep is None else sep.min(d)
                if not sep > 0:
                    lines = None
                    break
                sin = (res/sep).upper()
                theta = arb.pi()/2 if not sin < 1 else sin.asin().upper()
                lines[k, i] = (col, nv, theta)
            if lines is None:
                break
        out["lines" if tag == "box" else "plines"] = lines
        out["sel" if tag == "box" else "psel"] = select(s, H) if lines is not None else None
    return out


def select(s, H):
    """两步读数最小的三个是否为 {(1,j1),(2,j2),(0,c)}，c 取指定候选之一。返回各扇区的可能编号。"""
    j1, j2 = (1, LIGHT[s][1][0]), (2, LIGHT[s][2][0])
    c0 = [(0, j) for j in LIGHT[s][0]]
    S = {j1, j2, *c0}
    top = H[j1].upper().max(H[j2].upper())
    if len(c0) == 1:
        T = top.max(H[c0[0]].upper())
    else:
        # 至少一个扇区 0 候选严格重于 j1、j2，于是它是第四个；另一个候选是第三个。
        if not top < H[c0[0]].lower().max(H[c0[1]].lower()):
            return None
        T = top.max(H[c0[0]].upper().min(H[c0[1]].upper()))
    if not all(T < H[x].lower() for x in H if x not in S):
        return None
    if len(c0) == 1:
        zero = [c0[0][1]]
    else:
        low = H[c0[0]].upper().min(H[c0[1]].upper())
        zero = [j for _, j in c0 if not H[0, j].lower() > low]
    return {0: zero, 1: [j1[1]], 2: [j2[1]]}


DESIGNATED = {s: {k: LIGHT[s][k] for k in range(3)} for s in ABS}
CACHE = {}


def charge_data(q, x, y, h):
    s = abs(q)
    key = (s, (s*x) % 1, (s*y) % 1, s*h)
    if key not in CACHE:
        # 联合盒 |X - 2 pi x|, |Y - 2 pi y| <= 2 pi h；荷 s 读 e^{isX}，|dU| <= s * 2 pi h。
        CACHE[key] = analyze(s, turn((s*x) % 1), turn((s*y) % 1), radians(s*h))
    d = CACHE[key]
    return d if q > 0 else conjugate_data(d)


def conjugate_data(d):
    """L_{-s}(X,Y) = conj L_s(X,Y)：扇区 k 的数据来自扇区 -k，向量取共轭。"""
    out = {}
    for tag in ("lines", "plines"):
        v = d[tag]
        out[tag] = None if v is None else {((-k) % 3, i): ([z.conjugate() for z in col], nv, th)
                                           for (k, i), (col, nv, th) in v.items()}
    for tag in ("sel", "psel"):
        v = d[tag]
        out[tag] = None if v is None else {(-k) % 3: v[k] for k in range(3)}
    return out


def candidates(q, d, point):
    """盒上可能的轻线：选择已认证就用它；否则（仅盒）用阶段 A 的定理给出的全部指定候选。"""
    sel = d["psel" if point else "sel"]
    lines = d["plines" if point else "lines"]
    if lines is None:
        return None
    if sel is None:
        if point:
            return None
        sel = {k: LIGHT[abs(q)][k if q > 0 else (-k) % 3] for k in range(3)}
    return {k: [lines[k, i] for i in sel[k]] for k in range(3)}


def cos2_range(la, lb):
    lo, hi = None, None
    for ca, na, ta in la:
        for cb, nb, tb in lb:
            ip = sum((u.conjugate()*w for u, w in zip(ca, cb)), acb(0))
            cs = (abs(ip)/(na*nb)).min(arb(1))
            phi = cs.acos()
            pmin = (phi-ta-tb).lower().max(arb(0))
            pmax = (phi+ta+tb).upper().min(arb.pi()/2)
            vhi, vlo = (pmin.cos()**2).upper(), (pmax.cos()**2).lower()
            lo = vlo if lo is None else lo.min(vlo)
            hi = vhi if hi is None else hi.max(vhi)
    return lo, hi


def leakage_range(cq, cr):
    lo, hi, secs = arb(1), arb(1), []
    for k in range(3):
        cl, ch = cos2_range(cq[k], cr[k])
        lo -= ch/3
        hi -= cl/3
        secs.append((cl, ch))
    return lo.lower(), hi.upper(), secs


# ---------------- 阶段 A：每个荷在自身整张环面上的轻线标签 ----------------
def stage_a(s, g0=16, hmin=Fr(1, 2**12)):
    """四叉树覆盖荷 s 的自身环面，每个盒上认证两步最轻三个读数的扇区标签。"""
    stack = [(Fr(2*i+1, 2*g0), Fr(2*j+1, 2*g0), Fr(1, 2*g0)) for i in range(g0) for j in range(g0)]
    boxes = ambiguous = 0
    while stack:
        x, y, h = stack.pop()
        d = analyze(s, turn(x), turn(y), radians(h))
        if d["sel"] is not None:
            boxes += 1
            ambiguous += len(d["sel"][0]) > 1
        else:
            assert h > hmin
            h2 = h/2
            stack += [(x+a, y+b, h2) for a in (-h2, h2) for b in (-h2, h2)]
    return boxes, ambiguous


# ---------------- 阶段 B：两两泄漏在联合环面上的两侧界 ----------------
def stage_b(p, inf_t, sup_t, g0=16, hmin=Fr(1, 2**11)):
    """四叉树覆盖联合环面；每个盒上泄漏的区间须落在 [inf_t, sup_t] 内，否则四分。"""
    q, r = p
    inf_t, sup_t = arb(inf_t.numerator)/inf_t.denominator, arb(sup_t.numerator)/sup_t.denominator
    stack = [(Fr(2*i+1, 2*g0), Fr(2*j+1, 2*g0), Fr(1, 2*g0)) for i in range(g0) for j in range(g0)]
    st = dict(lo=None, hi=None, n=0, sec=[[None, None] for _ in range(3)])
    while stack:
        x, y, h = stack.pop()
        cq = candidates(q, charge_data(q, x, y, h), False)
        cr = candidates(r, charge_data(r, x, y, h), False)
        ok = cq is not None and cr is not None
        if ok:
            lo, hi, secs = leakage_range(cq, cr)
            ok = lo >= inf_t and hi <= sup_t
        if not ok:
            assert h > hmin, (p, x, y, h)
            h2 = h/2
            stack += [(x+a, y+b, h2) for a in (-h2, h2) for b in (-h2, h2)]
            continue
        st["n"] += 1
        st["lo"] = lo if st["lo"] is None else st["lo"].min(lo)
        st["hi"] = hi if st["hi"] is None else st["hi"].max(hi)
        for k in range(3):
            st["sec"][k][0] = secs[k][0] if st["sec"][k][0] is None else st["sec"][k][0].min(secs[k][0])
            st["sec"][k][1] = secs[k][1] if st["sec"][k][1] is None else st["sec"][k][1].max(secs[k][1])
    return st


def point_value(p, x, y):
    """联合相位 (2 pi x, 2 pi y) 处、对全部 delta 的泄漏区间；要求该点轻线选择已认证。"""
    q, r = p
    cq = candidates(q, charge_data(q, x, y, Fr(0)), True)
    cr = candidates(r, charge_data(r, x, y, Fr(0)), True)
    assert cq is not None and cr is not None and all(len(cq[k]) == len(cr[k]) == 1 for k in range(3))
    lo, hi, secs = leakage_range(cq, cr)
    return lo, hi, secs


def flagship_check():
    """在旗舰零点上与原定义（十九维全矩阵 Rump 本征向量、h < 9）逐对比较。"""
    light = run_path(str(ROOT / "ch39-local-yukawa-plucker.py"))
    old = ctx.prec
    ctx.prec = 96
    hgt = acb.zeta_zero(2_306_888_100_085).imag
    delta = arb(23)/25*arb(2)**-20/(hgt*hgt+SOFT["point"]["study"]["offset"])
    assert delta < DMAX
    bases = {q: light["light_basis"](q, hgt, 15-delta/3) for q in (1, 2, -4)}
    ph = GEO["phases"](hgt)
    out = {}
    for q, r in ((1, 2), (1, -4)):
        ov = dagger(bases[q])*bases[r]
        ref = 1-sum((abs(ov[i, j])**2 for i in range(3) for j in range(3)), arb(0))/3
        data = {}
        for c in (q, r):
            s = abs(c)
            d = analyze(s, ph[724]**s, ph[804]**s, arb(0))
            d = d if c > 0 else conjugate_data(d)
            data[c] = candidates(c, d, True)
            assert all(len(data[c][k]) == 1 for k in range(3))
        lo, hi, _ = leakage_range(data[q], data[r])
        assert lo <= ref.upper() and ref.lower() <= hi             # 两种算法给出同一个数
        out[q, r] = (ref, lo, hi)
    ctx.prec = old
    return out


def light_projector_h9(q, hgt, delta, x, y):
    """原定义：十九维全矩阵 Rump 本征向量，取 h < 9 的模，要求恰三个且无一个本征值的 h 区间跨过 9。"""
    lap = SHELL["laplacian"](q, hgt, phases={724: turn(x), 804: turn(y)})
    eig, vec = lap.eig(right=True, algorithm="rump")
    assert all(not z.overlaps(w) for j, z in enumerate(eig) for w in eig[:j])
    cols = []
    for j, z in enumerate(eig):
        h = abs(z.real*(z.real-(15-delta/3)))
        assert h < 9 or h > 9
        if h < 9:
            nrm = sum((abs(vec[i, j])**2 for i in range(19)), arb(0)).sqrt()
            cols.append([vec[i, j]/nrm for i in range(19)])
    assert len(cols) == 3
    B = acb_mat(list(zip(*cols)))
    return B*dagger(B)


def permutation_functional_upper(P, Q):
    """调相—重排下界泛函 1 - max_pi sum_i <sort|P_i|, sort|Q_{pi i}|>/3 的上界。

    取任一置换 pi 与任一行内配对都只会减小该和，故以 |元素| 的下端点、按数值排序后的配对
    与数值指派给出的 pi 求和，得到 max 的下界，即泛函的上界。"""
    from scipy.optimize import linear_sum_assignment
    rows = []
    for M in (P, Q):
        rows.append([sorted((abs(M[i, j]).lower() for j in range(19)), key=float, reverse=True)
                     for i in range(19)])
    sizes = [[sum((a*b for a, b in zip(u, v)), arb(0)) for v in rows[1]] for u in rows[0]]
    _, perm = linear_sum_assignment(-np.array([[float(z.mid()) for z in row] for row in sizes]))
    total = sum((sizes[i][perm[i]] for i in range(19)), arb(0))
    return (1-total/3).upper()


def permutation_functional_points():
    """旗舰零点的 delta 下，两处有理联合相位上 (1,2)、(1,-4) 的调相—重排下界泛函。"""
    old = ctx.prec
    ctx.prec = 96
    hgt = acb.zeta_zero(2_306_888_100_085).imag
    delta = arb(23)/25*arb(2)**-20/(hgt*hgt+SOFT["point"]["study"]["offset"])
    out = {}
    for (q, r), (x, y) in (((1, 2), (Fr(57, 64), Fr(35, 64))), ((1, -4), (Fr(33, 64), Fr(13, 64)))):
        out[q, r] = permutation_functional_upper(light_projector_h9(q, hgt, delta, x, y),
                                                 light_projector_h9(r, hgt, delta, x, y))
    ctx.prec = old
    return out


def handover_points():
    """偶荷扇区 0 的两个候选各在某点唯一被选中（自身相位，单位为圈；逐点模式，delta 取遍）。"""
    pts = {2: ((Fr(0), Fr(0), 0), (Fr(4, 64), Fr(30, 64), 3)),
           4: ((Fr(0), Fr(0), 0), (Fr(0), Fr(7, 16), 3)),
           6: ((Fr(0), Fr(0), 0), (Fr(0), Fr(3, 8), 3))}
    for s, cases in pts.items():
        for x, y, j in cases:
            d = analyze(s, turn(x), turn(y), arb(0))
            assert d["psel"] is not None and d["psel"][0] == [j], (s, x, y, d["psel"])
    return pts


# 每对：(证书下界门槛, 证书上界门槛, 下确界见证点, 上确界见证点)；见证点以圈为单位。
ODD_EVEN = [(1, 2), (1, -4), (1, 6), (2, -3), (-3, -4), (-3, 6)]
PLAN = {
    (1, 2): (Fr(2, 5), Fr(13, 20), (Fr(1, 2), Fr(1, 2)), (Fr(7505, 10**4), Fr(2502, 10**4))),
    (1, -4): (Fr(2, 5), Fr(13, 20), (Fr(1, 2), Fr(1, 2)), (Fr(8748, 10**4), Fr(1248, 10**4))),
    (1, 6): (Fr(2, 5), Fr(13, 20), (Fr(1, 2), Fr(1, 2)), (Fr(834, 10**4), Fr(5834, 10**4))),
    (2, -3): (Fr(2, 5), Fr(13, 20), (Fr(1, 2), Fr(1, 2)), (Fr(2659, 10**4), Fr(7419, 10**4))),
    (-3, -4): (Fr(2, 5), Fr(13, 20), (Fr(1, 2), Fr(1, 2)), (Fr(6279, 10**4), Fr(3754, 10**4))),
    (-3, 6): (Fr(2, 5), Fr(13, 20), (Fr(53640, 10**5), Fr(80065, 10**5)), (Fr(9205, 10**4), Fr(795, 10**4))),
    (1, -3): (Fr(1, 50), Fr(7, 100), (Fr(1, 2), Fr(0)), (Fr(1, 2), Fr(1, 2))),
    (2, -4): (Fr(1, 1000), Fr(9, 20), (Fr(5005, 10**4), Fr(3347, 10**4)), (Fr(3, 4), Fr(3, 4))),
    (2, 6): (Fr(1, 1000), Fr(9, 20), (Fr(33, 10**4), Fr(4964, 10**4)), (Fr(1471, 10**4), Fr(1474, 10**4))),
    (-4, 6): (Fr(1, 1000), Fr(9, 20), (Fr(4822, 10**4), Fr(5178, 10**4)), (Fr(1509, 10**4), Fr(6511, 10**4))),
}


# 证书数值（见说明“关键数字”）：盒证书的下确界下界、上确界上界；见证点上的泄漏上界与下界。
INF_LO = {(1, 2): ".40006", (1, -4): ".40011", (1, 6): ".43861", (2, -3): ".40045", (-3, -4): ".400007",
          (-3, 6): ".42822", (1, -3): ".023740", (2, -4): ".0010167", (2, 6): ".0010002", (-4, 6): ".0027328"}
SUP_HI = {(1, 2): ".60810", (1, -4): ".60772", (1, 6): ".58294", (2, -3): ".64999", (-3, -4): ".64039",
          (-3, 6): ".60428", (1, -3): ".069975", (2, -4): ".39162", (2, 6): ".42499", (-4, 6): ".40179"}
WIT_MIN = {(1, 2): ".5166743", (1, -4): ".5242890", (1, 6): ".5231723", (2, -3): ".5219934",
           (-3, -4): ".5128896", (-3, 6): ".5324159", (1, -3): ".04262503", (2, -4): ".008586716",
           (2, 6): ".003338585", (-4, 6): ".02123601"}
WIT_MAX = {(1, 2): ".5345892", (1, -4): ".5383288", (1, 6): ".5353442", (2, -3): ".5422583",
           (-3, -4): ".5356564", (-3, 6): ".5487335", (1, -3): ".04693149", (2, -4): ".3493421",
           (2, 6): ".3455486", (-4, 6): ".3612111"}

STAGE_A = {1: (2014, 0), 2: (940, 702), 3: (2512, 0), 4: (988, 700), 6: (1462, 1050)}
STAGE_B = {(1, 2): 3346, (1, -4): 12898, (1, 6): 16384, (2, -3): 4180, (-3, -4): 14848, (-3, 6): 16384,
           (1, -3): 7870, (2, -4): 15400, (2, 6): 68392, (-4, 6): 16384}


def verify():
    t0 = time.process_time()
    flag = flagship_check()
    for p, (ref, lo, hi) in flag.items():
        print("旗舰零点", p, "原定义", ref.str(10), "扇区分解", lo.str(10), hi.str(10))
    assert flag[1, 2][0] > arb(".51702520") and flag[1, -4][0] > arb(".52481329")
    assert flag[1, 2][0] < arb(".51702521") and flag[1, -4][0] < arb(".52481330")
    assert flag[1, 2][0]-arb(WIT_MIN[1, 2]) > arb("3.5e-4")       # 旗舰值高于 X = Y = pi 处的值
    pf = permutation_functional_points()
    print("调相—重排下界泛函的上界", {p: v.str(12) for p, v in pf.items()})
    assert pf[1, 2] < arb(".30145") and pf[1, -4] < arb(".30553")
    handover_points()
    for s in ABS:
        boxes, amb = stage_a(s)
        print("阶段 A 荷", s, "盒数", boxes, "扇区零候选未分开的盒", amb)
        assert (boxes, amb) == STAGE_A[s]
    results = {}
    for p, (inf_t, sup_t, wmin, wmax) in PLAN.items():
        st = stage_b(p, inf_t, sup_t)
        a_lo, a_hi, _ = point_value(p, *wmin)
        b_lo, b_hi, _ = point_value(p, *wmax)
        results[p] = (st, a_hi, b_lo)
        print("阶段 B", p, "盒数", st["n"], "下确界 ∈", st["lo"].str(6), a_hi.str(8),
              "上确界 ∈", b_lo.str(8), st["hi"].str(6),
              "扇区 cos^2 范围", [(x.str(4), y.str(4)) for x, y in st["sec"]])
        assert st["lo"] >= arb(inf_t.numerator)/inf_t.denominator
        assert st["hi"] <= arb(sup_t.numerator)/sup_t.denominator and a_hi < b_lo
        assert st["lo"] >= arb(INF_LO[p]) and st["hi"] < arb(SUP_HI[p])
        assert a_hi < arb(WIT_MIN[p]) and b_lo > arb(WIT_MAX[p]) and st["n"] == STAGE_B[p]
        if p == (1, 2):
            assert a_lo > arb(".516674")
        if p in ODD_EVEN:                     # 盒区间之并含全部真值：扇区 0 近乎正交，扇区 1、2 部分重合
            assert st["sec"][0][1] < arb(".22")
            assert all(st["sec"][k][0] > arb(".37") and st["sec"][k][1] < arb(".98") for k in (1, 2))
    print("用时", round(time.process_time()-t0, 1), "秒")
    return flag, results


if __name__ == "__main__":
    verify()
