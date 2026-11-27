"""二次增长裕度选权在一族弦度量上的度量依赖：整族严格真空、w_105 与 w_104 的包络、相对效率极大极小权。

模型对象沿用 ch39-quadratic-growth-margin-vacuum.py：105 圈切片 P（旗舰相位 t0=0.01046，A=(x,y)），
N_w(B)=S_w(B;0)-S_w(A;0)。对度量系数 c=(c1,c2,c3)（c1,c2>0，c3>=0）取弦距离平方
    rho_c(B)^2 = min_{C=±A} [c1 r1_C(B) + c2 r2_C(B) + c3 r3_C(B)],  r_C = 2-2cos k.(B-C)，
k 依次为 (1,0)、(0,1)、(1,-1)；裕度 g_c(w)=inf_{B≠±A} N_w(B)/rho_c(B)^2，G(c)=max_P g_c。
两种已用度量是 c=(1,1,1) 与 c=(1,1,0)；它们张成的锥按比例归一为 c_s=(1,1,s)，s in [0,1]。

覆盖在已有的二阶泰勒格上增加一类折点叶格：h=max_C h_C >= lam h_A+(1-lam) h_{-A}。
本脚本认证：
  1. 漂移证书：凡 g_{(1,1,0)}(w) >= 1/1000 的 w in P，镜像漂移 < m0 < 0；代入全高度引理，这些 w
     都是全部高度上的严格共同真空（由此推出整族度量的最优权都是严格真空）。
  2. 节点 s=k/20：有理见证 W_k 与全环面覆盖给出 G(s_k) >= L_k，对偶证书给出 G(s_k) <= U_k，
     并包住 M(s_k) 的 w_105、w_104。
  3. 单元 [s_k, s_{k+1}]：由 G 单调与 sG 单调，M(s) 的每个元满足节点行，包住 w_105、w_104。
  4. 相对效率 e(w)=min_s g_s(w)/G(s) 的两侧界由 ch39-margin-metric-maximin.py 认证，证书在本 .json 的
     "maximin" 段。
浮点只用于挑选见证、竞争点与乘子；全部不等式由 Arb 球算术认证（精度 160 位）。
"""
from functools import cache
from pathlib import Path
from runpy import run_path
import json
import sys

from flint import acb, arb

PREC = 160


def cover(hessian, modes, x, y, gam, sig, box=arb("0.00002"), limit=2**24):
    """证明 h(B)=N(B)-min_C[gam(r1_C+r2_C)+sig r3_C] >= 0 于全环面（gam, sig >= 0）。

    h >= h_C := N - gam(r1_C+r2_C) - sig r3_C，h_C 沿三条轴可分。每条轴的一维函数取格中心的值与导数、
    格上二阶导数的球界，合成二维二阶泰勒下界。叶格四类：某个 h_C 的下界为正；对按中心梯度选取的有理
    lam in [0,1]，lam h_A+(1-lam) h_{-A} 的下界为正（h=max_C h_C 不小于它）；落在 ±A 的方框内；完全落在 rho_C < rho_{-C} 一侧且 grad h_C 不含零（此处 h=h_C，不可能是 h 的极小点）。
    方框内 Hess h_C >= Hess N - 2[[gam+sig,-sig],[-sig,gam+sig]] 正定，h_C(C)=0、grad h_C(C)=0。"""
    centre = (x, y, x-y)
    wt = (gam, gam, sig)
    trial = arb(0)
    for q, cs in modes:
        for c, t in zip(cs, centre):
            trial += (c*(1-acb(0, q*t).exp())).real

    @cache
    def axis_data(axis, n, index):
        if axis == 2:
            tc, r = arb.pi()*arb(2*index)/n, 2*arb.pi()/n
        else:
            tc, r = arb.pi()*(arb(2*index+1)/n-1), arb.pi()/n
        t = tc+arb(0, r)
        f = f1 = f2 = arb(0)
        for q, cs in modes:
            c = cs[axis]
            rc, ri = c*acb(0, q*tc).exp(), c*acb(0, q*t).exp()
            f += c.real-rc.real
            f1 += q*rc.imag
            f2 += q*q*ri.real
        out = []
        for sg in (1, -1):
            a = sg*centre[axis]
            uc, ui = tc-a, t-a
            lw = wt[axis]
            g, g1, g2 = f-lw*(2-2*uc.cos()), f1-2*lw*uc.sin(), f2-2*lw*ui.cos()
            out.append((g, g1, arb(min(0, g2.lower()))*r*r/2, g1+arb(0, (abs(g2)*r).upper()), lw*(2-2*ui.cos()), g2, r))
        return out

    def lower_mix(hs, n):
        # h = max_C h_C >= lam h_A + (1-lam) h_{-A}（lam in [0,1]）；lam 取使中心合成梯度最小的有理数。
        r = arb.pi()/n
        ga = [float((hs[0][0][1]+hs[2][0][1]).mid()), float((hs[1][0][1]-hs[2][0][1]).mid())]
        gm = [float((hs[0][1][1]+hs[2][1][1]).mid()), float((hs[1][1][1]-hs[2][1][1]).mid())]
        dx, dy = ga[0]-gm[0], ga[1]-gm[1]
        den = dx*dx+dy*dy
        lamf = 0.5 if den == 0 else min(1.0, max(0.0, -(gm[0]*dx+gm[1]*dy)/den))
        lam = arb(round(lamf*2**20))/2**20
        mu = 1-lam
        val = sum(lam*h[0][0]+mu*h[1][0]+arb(min(0, (lam*h[0][5]+mu*h[1][5]).lower()))*h[0][6]*h[0][6]/2 for h in hs)
        gx = lam*(hs[0][0][1]+hs[2][0][1])+mu*(hs[0][1][1]+hs[2][1][1])
        gy = lam*(hs[1][0][1]-hs[2][0][1])+mu*(hs[1][1][1]-hs[2][1][1])
        return val-(abs(gx)+abs(gy))*r

    def lower(hs, s, n):
        r = arb.pi()/n
        gx, gy = hs[0][s][1]+hs[2][s][1], hs[1][s][1]-hs[2][s][1]
        return sum(h[s][0]+h[s][2] for h in hs)-(abs(gx)+abs(gy))*r

    pending, counts = [(1, 0, 0)], [0, 0, 0, 0]
    while pending:
        n, i, j = pending.pop()
        hs = (axis_data(0, n, i), axis_data(1, n, j), axis_data(2, n, i-j))
        # 格中心处 h = max_C h_C；两者都为负即找到 h < 0 的点，证书不成立。
        assert not all(sum(h[s][0] for h in hs)-trial < 0 for s in (0, 1)), ("negative", n, i, j)
        if any(lower(hs, s, n)-trial > 0 for s in (0, 1)):
            counts[0] += 1
            continue
        if lower_mix(hs, n)-trial > 0:
            counts[3] += 1
            continue
        a, b = [arb.pi()*(arb(2*k+1, 1)/n-1) for k in (i, j)]
        if any(abs(a-s*x) < box and abs(b-s*y) < box for s in (-1, 1)):
            counts[1] += 1
            continue
        ra, rm = sum(h[0][4] for h in hs), sum(h[1][4] for h in hs)
        excluded = False
        for s, side in ((0, ra < rm), (1, rm < ra)):
            if side:
                px, py = hs[0][s][3]+hs[2][s][3], hs[1][s][3]-hs[2][s][3]
                excluded |= not px.contains(0) or not py.contains(0)
        if excluded:
            counts[2] += 1
            continue
        assert n < limit, (n, i, j, a, b)
        pending.extend((2*n, 2*i+u, 2*j+v) for u in (0, 1) for v in (0, 1))
    for s in (-1, 1):
        hxx, hxy, hyy = hessian(modes, s*x+arb(0, box), s*y+arb(0, box))
        p, r, d = hxx-2*(gam+sig), hyy-2*(gam+sig), hxy+2*sig
        assert p > 0 and p*r > d*d
    return counts


class Rows:
    """竞争点行与 A 处方向行：对一切 w in P 与一切 s >= 0，只要 g_s(w) >= 0，
    sum_q w_q a_q >= g_s(w)·min_C(R12_C + s R3_C)。
    竞争点行：a_q = S_q(B;0)-S_q(A;0)，对 C=±A 记 (R12_C, R3_C) = (r1_C+r2_C, r3_C)；
    方向行：a_q = u'Hess S_q(A;0)u/2，R12 = |u|^2，R3 = (u1-u2)^2（方向极限，两个 C 相同）。
    若 0 <= alpha <= g 且 0 <= beta <= s g，则 g·min_C(R12_C + s R3_C) >= min_C(alpha R12_C + beta R3_C)。"""

    def __init__(self, model, spec):
        self.model = model
        self.data = []
        for kind, u, v in spec:
            u, v = arb(u), arb(v)
            if kind == "p":
                a = {q: model.lib["energy"](q, model.bare[q], u, v)-model.home[q] for q in range(1, 106)}
                r = []
                for sg in (1, -1):
                    d1, d2 = u-sg*model.x, v-sg*model.y
                    r.append((2-2*d1.cos()+2-2*d2.cos(), 2-2*(d1-d2).cos()))
                    assert r[-1][0] > 0
            else:
                a = {}
                for q in range(1, 106):
                    _, _, hxx, hxy, hyy = model.lib["derivatives"]([(q, model.bare[q])], model.x, model.y)
                    a[q] = (u*u*hxx+2*u*v*hxy+v*v*hyy)/2
                r = [(u*u+v*v, (u-v)**2)]*2
                assert r[0][0] > 0
            self.data.append((a, r))

    def level(self, pairs):
        """pairs: [(行号, alpha, beta)]，给出 (a_q, b=min_C(alpha R12_C + beta R3_C))。"""
        out = []
        for j, al, be in pairs:
            a, r = self.data[j]
            out.append((a, (al*r[0][0]+be*r[0][1]).min(al*r[1][0]+be*r[1][1])))
        return out


def frac(text):
    p, q = text.split("/")
    return arb(int(p))/int(q)


def verify(root):
    base = run_path(str(root/"ch39-quadratic-growth-margin-vacuum.py"))
    model = base["Model"](root)
    hessian = base["hessian"]
    data = json.loads(Path(__file__).with_suffix(".json").read_text())
    model.set_phase(arb(data["t"]))
    rows = Rows(model, data["rows"])
    unit = {105: arb(1)}, {105: arb(-1)}, {104: arb(1)}, {104: arb(-1)}

    def enclose(certs, pairs):
        lv = rows.level(pairs)
        out = []
        for key, phi in zip(("w105", "-w105", "w104", "-w104"), unit):
            out.append(model.dual_value(certs[key], lv, phi, arb(1))[0])
        return (-out[1], out[0]), (-out[3], out[2])

    # ---- 1. 漂移证书与全高度引理：g_{(1,1,0)}(w) >= l0 推出严格共同真空。
    l0 = arb(1)/1000
    pairs0 = [(j, l0, arb(0)) for j in range(len(rows.data))]
    drift0 = model.dual_value(data["drift0"], rows.level(pairs0), model.mu, arb(1))[0]
    m0 = arb("-1.46")
    assert drift0 < m0
    hs, G = base["strict_at_all_heights"](model, l0, m0, 2)
    assert hs < arb("66151") and G < arb("362.24")
    print("漂移证书：g_(1,1,0) >= 1/1000 的权的镜像漂移上界、全高度引理常数 H_s、G", drift0.str(6), hs.str(6), G.str(6), flush=True)

    # ---- 2. 节点。
    nodes = data["nodes"]
    N = len(nodes)-1
    assert N == 20 and [nd["s"] for nd in nodes] == [f"{k}/20" for k in range(21)]
    S, L, U, E105, E104, supports = [], [], [], [], [], []
    for k, nd in enumerate(nodes):
        s = frac(nd["s"])
        W = model.complete({int(q): arb(v) for q, v in nd["free"].items()}, nd["basis"])
        modes = [(q, [W[q]*c for c in model.bare[q]]) for q in W]
        Lk = arb(nd["L"])
        cover(hessian, modes, model.x, model.y, Lk, Lk*s)
        allpairs = [(j, arb(1), s) for j in range(len(rows.data))]
        numer, norm = model.dual_value(nd["upper"], rows.level(allpairs), {}, arb(0))
        assert norm > 0
        Uk = numer/norm
        assert Lk < Uk
        e105, e104 = enclose(nd["enc"], [(j, Lk, Lk*s) for j in range(len(rows.data))])
        S.append(s); L.append(Lk); U.append(Uk); E105.append(e105); E104.append(e104)
        supports.append(set(W))
        print("节点 s、L <= G(s) <= U、支撑、M(s) 的 w_105 与 w_104 区间", nd["s"], Lk.str(8), Uk.str(10), len(W),
              e105[0].str(7), e105[1].str(7), e104[0].str(7), e104[1].str(7), flush=True)
    # 节点见证的支撑：11 到 14 项，共同部分恰为九圈。
    assert all(11 <= len(x) <= 14 for x in supports)
    assert set.intersection(*supports) == {1, 32, 35, 38, 40, 67, 81, 104, 105}
    # 与已有证书相容：G(1) > 3.0844、G(0) > 6.3487。
    assert L[-1] > arb("3.0844") and L[0] > arb("6.3487")
    # 推论的比值界：c 的最大分量不超过 min(c1,c2) 的 3000 倍时 min(c1,c2) G(c) >= L(1)/3000 > 1/1000。
    assert L[-1]/3000 > l0
    # G 单调不增、sG 单调不减在节点读数上不被违反（独立一致性检查）。
    for k in range(N):
        assert L[k+1] < U[k] and S[k]*L[k] < S[k+1]*U[k+1]

    # ---- 3. 单元：s in [s_k, s_{k+1}]，w in M(s) 满足 a_j w >= L_{k+1} R12 + s_k L_k R3。
    C105, C104 = [], []
    for k in range(N):
        cert = data["cells"][k]
        e105, e104 = enclose(cert, [(j, L[k+1], S[k]*L[k]) for j in range(len(rows.data))])
        C105.append(e105)
        C104.append(e104)

    def reduce(vals, op):
        best = vals[0]
        for v in vals[1:]:
            best = op(best, v)
        return best
    out105 = (reduce([c[0].lower() for c in C105+E105], lambda a, b: a.min(b)),
              reduce([c[1].upper() for c in C105+E105], lambda a, b: a.max(b)))
    out104 = (reduce([c[0].lower() for c in C104+E104], lambda a, b: a.min(b)),
              reduce([c[1].upper() for c in C104+E104], lambda a, b: a.max(b)))
    in105 = (reduce([e[1].upper() for e in E105], lambda a, b: a.min(b)),
             reduce([e[0].lower() for e in E105], lambda a, b: a.max(b)))
    in104 = (reduce([e[1].upper() for e in E104], lambda a, b: a.min(b)),
             reduce([e[0].lower() for e in E104], lambda a, b: a.max(b)))
    # 数字断言：外包络、内包络（内包络由连通性给出区间），以及选取的非单调性。
    assert out105[0] > arb("0.10474") and out105[1] < arb("0.12008")
    assert out104[0] > arb("0.022296") and out104[1] < arb("0.029174")
    assert in105[0] < arb("0.110772") and in105[1] > arb("0.112738")
    assert in104[0] < arb("0.023395") and in104[1] > arb("0.026046")
    table = {0: ("6.348723", "6.3487299", "0.110870", "0.110874", "0.0247775", "0.0247784"),
             2: ("5.680244", "5.6802501", "0.112367", "0.112381", "0.0247533", "0.0247810"),
             5: ("5.022098", "5.0221987", "0.110690", "0.110772", "0.0252885", "0.0253450"),
             11: ("3.912629", "3.9127081", "0.112738", "0.112967", "0.0253894", "0.0255828"),
             20: ("3.084429", "3.0844914", "0.111185", "0.111366", "0.0232412", "0.0233950")}
    for k, (l, u, a, b, c, d) in table.items():
        assert nodes[k]["L"] == l and U[k] < arb(u)
        assert arb(a) < E105[k][0] and E105[k][1] < arb(b) and arb(c) < E104[k][0] and E104[k][1] < arb(d)
    # 全部节点最优集的 w_105 上端都小于 0.112967（相对效率见证的 w_105 在其上方，见 ch39-margin-metric-maximin.py）。
    assert all(e[1] < arb("0.112967") for e in E105)
    # 任取 w(s) in M(s)：w_105(0) < w_105(1/10) > w_105(1/4) < w_105(11/20) > w_105(1)。
    chain = [0, 2, 5, 11, 20]
    for p, q, up in zip(chain, chain[1:], (True, False, True, False)):
        assert (E105[p][1] < E105[q][0]) if up else (E105[q][1] < E105[p][0])
    print("w_105 外包络、内包络", out105[0].str(7), out105[1].str(7), in105[0].str(7), in105[1].str(7))
    print("w_104 外包络、内包络", out104[0].str(7), out104[1].str(7), in104[0].str(7), in104[1].str(7), flush=True)
    return dict(out105=out105, out104=out104, in105=in105, in104=in104, L=L, U=U, rows=rows, model=model,
                data=data, S=S, hessian=hessian, E105=E105, E104=E104, base=base)


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else here)
