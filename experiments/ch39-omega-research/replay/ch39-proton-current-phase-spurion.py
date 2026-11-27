"""辅助电流重叠 c=g_L^dag g_D 的相位：Q-L 接触探针中为纯约定，加入 Delta B=3 单点顶点后成为可观测量。
配套说明：辅助电流重叠的相位在接触探针中可吸收而在重子数改变3的单点替身顶点下进入归一无关读数.md

(1) 探针类中的冗余：Q -> e^{-i theta/3} Q（共轭生成元取相反相位）是格拉斯曼代数的自同构，Berezin 雅可比为 1，
    把 K21=Fc 变成 |Fc|；带 U(1)_Q 荷 q 的算符乘 e^{-i q theta/3}。脚本核对各双线性与各耦合的荷。
(2) 完整模型的重相格（书面证明的整数线性代数部分）：在汤川矩阵有代间混合与弱瞬子顶点下，
    不变量 n arg(Fc) + m arg(eta_I) + (味相位) 须满足 n+3m=0；脚本给出显式整数组合、重子数配对，
    并核对只有对角汤川时无含 arg(Fc) 的不变量、去掉 Weinberg 项时结论不变。
(3) 单点扩充探针（新的模型选择，单点替身顶点）：在接触项之外加 eta W^3 + conj(eta) Wbar^3，W=sum_a (QL)_a D_a，
    Wbar=sum_a Dbar_a T_a。Z 的精确多项式、读数比 X_AB, X_AA, X_BB 的闭式，以及 X_AB 随
    Phi = 3 arg(Fc) - arg(eta) 改变的精确有理证书。
全部为精确整数/有理运算。
"""
import runpy
import time
from collections import defaultdict
from pathlib import Path

import sympy as sp

HERE = Path(__file__).resolve().parent
_static = runpy.run_path(str(HERE / "ch39-proton-current-overlap-static.py"))


# ---------------------------------------------------------------- (1) U(1)_Q 与 U(1)_L 荷
def charges():
    """每个 Q 分量荷 (1,0)，L 分量荷 (0,1)；共轭生成元取相反荷。
    双线性的荷由其单项的生成元直接数出，接触项各耦合的荷为其算符荷的相反数。"""
    D, QL = _static["bilinears"](False)
    Db, T = _static["bilinears"](True)

    def charge(poly, barred):
        out = set()
        for m in poly:
            nq = sum(1 for k in range(12) if m >> k & 1)
            nl = sum(1 for k in range(12, 16) if m >> k & 1)
            out.add((-nq, -nl) if barred else (nq, nl))
        assert len(out) == 1
        return out.pop()
    ch = {}
    for name, polys, barred in [("D", D, False), ("QL", QL, False), ("Db", Db, True), ("T", T, True)]:
        cs = {charge(p, barred) for p in polys}
        assert len(cs) == 1
        ch[name] = cs.pop()
    assert ch == {"D": (2, 0), "QL": (1, 1), "Db": (-2, 0), "T": (-1, -1)}
    add = lambda a, b: (a[0]+b[0], a[1]+b[1])
    op = {"K11": add(ch["Db"], ch["D"]), "K22": add(ch["QL"], ch["T"]),
          "K21": add(ch["QL"], ch["D"]), "K12": add(ch["Db"], ch["T"])}
    assert op == {"K11": (0, 0), "K22": (0, 0), "K21": (3, 1), "K12": (-3, -1)}
    W3 = (3*op["K21"][0], 3*op["K21"][1])
    assert W3 == (9, 3)
    # 相位不变量 n arg K21 + m arg eta 须满足 n(3,1)+m(9,3)=0  <=>  n=-3m：唯一的生成元为 arg eta - 3 arg K21。
    n, m = sp.symbols("n m")
    sol = sp.solve([3*n+9*m, n+3*m], [n], dict=True)
    assert sol == [{n: -3*m}]
    return op


# ---------------------------------------------------------------- (2) 完整模型的重相格
def full_model_lattice():
    """场：Q_g,u_g,d_g,L_g,e_g (g=1,2,3)，H。只用算符荷向量（耦合荷为其相反数，不影响张成）。
    汤川 Q_g H u_h, Q_g H^dag d_h, L_g H^dag e_h，Weinberg L_g L_h H H，弱瞬子 prod_g Q_g^3 L_g，接触项 Q_1^3 L_1。
    相位组合 n arg(Fc) + m arg(eta_I) + (汤川、Weinberg 相位) 是重相不变量
      <=>  n v_c + m v_I 落在汤川与 Weinberg 向量的张成中
      <=>  对零化子中每个 q 有 q.(n v_c + m v_I) = 0。
    三种耦合布局：(a) 汤川全部 9 项非零 + Weinberg 全部 6 项非零；(b) 只有对角汤川 + 对角 Weinberg；
    (c) 汤川全部非零、无 Weinberg 项。"""
    names = [f"{f}{g}" for f in "QudLe" for g in (1, 2, 3)]+["H"]
    idx = {s: i for i, s in enumerate(names)}

    def vec(**kw):
        v = [0]*len(names)
        for k, c in kw.items():
            v[idx[k]] += c
        return sp.Matrix(v)
    yu = lambda g, h: vec(**{f"Q{g}": 1, "H": 1, f"u{h}": 1})
    yd = lambda g, h: vec(**{f"Q{g}": 1, "H": -1, f"d{h}": 1})
    ye = lambda g, h: vec(**{f"L{g}": 1, "H": -1, f"e{h}": 1})

    def wv(g, h):
        v = vec(**{f"L{g}": 1, "H": 2})
        v[idx[f"L{h}"]] += 1
        return v
    pairs = [(g, h) for g in (1, 2, 3) for h in (1, 2, 3)]
    yuk_full = [f(g, h) for g, h in pairs for f in (yu, yd, ye)]
    wein_full = [wv(g, h) for g in (1, 2, 3) for h in range(g, 4)]
    layouts = {
        "full": yuk_full+wein_full,
        "diagonal": [f(g, g) for g in (1, 2, 3) for f in (yu, yd, ye)]+[wv(g, g) for g in (1, 2, 3)],
        "no-weinberg": yuk_full,
    }
    inst = sum((vec(**{f"Q{g}": 3, f"L{g}": 1}) for g in (1, 2, 3)), sp.zeros(len(names), 1))
    cont = vec(Q1=3, L1=1)
    B = sp.Matrix([sp.Rational(1, 3) if s[0] == "Q" else (-sp.Rational(1, 3) if s[0] in "ud" else 0)
                   for s in names])
    Lnum = sp.Matrix([{"L": 1, "e": -1}.get(s[0], 0) for s in names])
    Yh = sp.Matrix([{"Q": sp.Rational(1, 6), "u": -sp.Rational(2, 3), "d": sp.Rational(1, 3), "L": -sp.Rational(1, 2),
                     "e": 1, "H": sp.Rational(1, 2)}[s[0]] for s in names])
    n, m = sp.symbols("n m")
    sols = {}
    for name, vecs in layouts.items():
        Y = sp.Matrix.hstack(*vecs)
        null = Y.T.nullspace()
        # (n,m) 的解空间
        eqs = [(q.T*(n*cont+m*inst))[0] for q in null]
        M = sp.Matrix([[sp.diff(e, n), sp.diff(e, m)] for e in eqs])
        sols[name] = (len(null), M.nullspace())
        assert (Y.T*B).is_zero_matrix and (Y.T*Yh).is_zero_matrix
    # (a) 零化子二维，由 B 与超荷张成；解为 n+3m=0。
    dim, ns = sols["full"]
    assert dim == 2 and len(ns) == 1 and ns[0][0] == -3*ns[0][1]
    # (b) 只有对角汤川：零化子四维（代别重子数），一切解都有 n=0，不存在含 arg(Fc) 的不变量。
    dim, ns = sols["diagonal"]
    assert dim == 4 and all(v[0] == 0 for v in ns)
    # (c) 去掉 Weinberg：零化子三维（B、轻子数、超荷），解仍是唯一方向 n+3m=0。
    dim, ns = sols["no-weinberg"]
    Yc = sp.Matrix.hstack(*layouts["no-weinberg"])
    assert dim == 3 and (Yc.T*Lnum).is_zero_matrix and sp.Matrix.hstack(*Yc.T.nullspace(), B, Lnum, Yh).rank() == 3
    assert len(ns) == 1 and ns[0][0] == -3*ns[0][1]
    # 重子数配对：B.v_c = 1, B.v_I = 3，汤川与 Weinberg 向量与 B 正交；整数组合 n v_c + m v_I + sum a_i y_i = 0
    # 配 B 得 n + 3m = 0（只用到 B 是这些项的对称）。
    assert (B.T*inst)[0] == 3 and (B.T*cont)[0] == 1 and (Yh.T*inst)[0] == 0 and (Yh.T*cont)[0] == 0
    # 显式整数实现 (n,m)=(3,-1)：3 v_c - v_I = sum_{g=2,3} [3(y^u_{11}-y^u_{g1}) + (w_{11}-w_{1g})]；
    # 无 Weinberg 时以 y^e 代替：3 v_c - v_I = sum_{g=2,3} [3(y^u_{11}-y^u_{g1}) + (y^e_{11}-y^e_{g1})]。
    combo = sum((3*(yu(1, 1)-yu(g, 1))+(wv(1, 1)-wv(1, g)) for g in (2, 3)), sp.zeros(len(names), 1))
    assert combo == 3*cont-inst
    combo_e = sum((3*(yu(1, 1)-yu(g, 1))+(ye(1, 1)-ye(g, 1)) for g in (2, 3)), sp.zeros(len(names), 1))
    assert combo_e == 3*cont-inst
    return names


# ---------------------------------------------------------------- (3) 单点扩充探针
def extended_site():
    D, QL = _static["bilinears"](False)
    Db, T = _static["bilinears"](True)
    tu, tb = _static["top_table"](D, QL), _static["top_table"](Db, T)
    k11, k22, k21, k12, eta, etab = sp.symbols("k11 k22 k21 k12 eta etab")
    # 偶元互相对易：把 D_a,(QL)_a,Dbar_a,T_a 当作交换符号展开，顶阶系数由 top_table（直接格拉斯曼枚举）给出，
    # 非共轭一侧与共轭一侧分解。符号 0..2 D_a，3..5 (QL)_a，6..8 Dbar_a，9..11 T_a。
    ok = lambda e: sum(e[0:3]) <= 4 and sum(e[3:6]) <= 4 and sum(e[6:9]) <= 4 and sum(e[9:12]) <= 4

    def pmul(A, B):
        out = defaultdict(lambda: 0)
        for ea, ca in A.items():
            for eb, cb in B.items():
                e = tuple(x+y for x, y in zip(ea, eb))
                if ok(e):
                    out[e] += ca*cb
        return {k: sp.expand(v) for k, v in out.items() if v != 0}

    def mono(*ix):
        e = [0]*12
        for i in ix:
            e[i] += 1
        return tuple(e)
    one = {tuple([0]*12): sp.Integer(1)}
    expE = one
    for a in range(3):
        for m, c in [(mono(6+a, a), k11), (mono(3+a, 9+a), k22), (mono(3+a, a), k21), (mono(6+a, 9+a), k12)]:
            f, cur = dict(one), one
            for nn in range(1, 5):
                cur = pmul(cur, {m: c})
                for k, v in cur.items():
                    f[k] = f.get(k, 0)+v/sp.factorial(nn)
            expE = pmul(expE, f)
    W = {mono(3+a, a): 1 for a in range(3)}
    Wb = {mono(6+a, 9+a): 1 for a in range(3)}
    W3, Wb3 = pmul(pmul(W, W), W), pmul(pmul(Wb, Wb), Wb)
    assert W3 and Wb3
    # (W^3)^2 与 W^3 Wbar^0 以外的高次项超出单侧次数 8，exp(eta W^3)=1+eta W^3 精确。
    assert pmul(W3, W3) == {}
    tot = pmul(pmul(expE, {**one, **{k: eta*v for k, v in W3.items()}}), {**one, **{k: etab*v for k, v in Wb3.items()}})
    Z = 0
    for e, c in tot.items():
        d, l, eb, m = e[0:3], e[3:6], e[6:9], e[9:12]
        if sum(d) == 4 and sum(l) == 4 and sum(eb) == 4 and sum(m) == 4:
            Z += c*tu.get((d, l), 0)*tb.get((eb, m), 0)
    Z = sp.expand(Z)
    rho = sp.Symbol("rho")
    p = 1+24*rho+90*rho**2+80*rho**3+15*rho**4
    s = k11*k22
    closed = 2211840*(s**4*p.subs(rho, k12*k21/s)+120*(eta*k12**3+etab*k21**3)*(4*s+3*k12*k21)
                      + 2880*eta*etab*(s+3*k12*k21))
    assert sp.expand(closed-Z) == 0
    # U(1)_Q 重相不变性：K21->z^3 K21, K12->K12/z^3, eta->z^9 eta, conj(eta)->conj(eta)/z^9。
    zz = sp.Symbol("zz")
    assert sp.expand(Z.subs({k21: zz**3*k21, k12: k12/zz**3, eta: zz**9*eta, etab: etab/zz**9}, simultaneous=True)-Z) == 0
    # 读数比：A、B 分别是 K11、K22 的算符。
    XAB = sp.diff(Z, k11, k22)*Z/(sp.diff(Z, k11)*sp.diff(Z, k22))
    XAA = sp.diff(Z, k11, 2)*Z/sp.diff(Z, k11)**2
    XBB = sp.diff(Z, k22, 2)*Z/sp.diff(Z, k22)**2
    # 不变量（K11,K22 取正实数，K12=conj(K21)，eta 的共轭为 etab）：rho=K12K21/(K11K22)，
    # mu=Re(eta K12^3)/(K11K22)^3，sigma=|eta|^2/(K11K22)^3。
    # 取 K11=K22=1，K12=K21=a=sqrt(rho)（Fc 取正实数，相位全在 eta），eta=h e^{i psi}，Phi=-psi。
    a, mu, sg = sp.symbols("a mu sigma", real=True)
    R = sp.Symbol("R")
    zt = p.subs(rho, R)+240*(4+3*R)*mu+2880*(1+3*R)*sg
    XAB_closed = (2+27*R+45*R**2+10*R**3+120*mu+360*sg)*zt/(2*(1+18*R+45*R**2+20*R**3+240*mu+720*sg)**2)
    hh, cc = sp.symbols("hh cc", real=True)  # eta = hh*(cc + i ss), ss^2 = 1-cc^2
    ss = sp.Symbol("ss", real=True)
    sub = {k11: 1, k22: 1, k12: a, k21: a, eta: hh*(cc+sp.I*ss), etab: hh*(cc-sp.I*ss)}
    back = {R: a**2, mu: hh*cc*a**3, sg: hh**2}
    num = sp.cancel(XAB.subs(sub))
    assert sp.simplify((num-XAB_closed.subs(back)).subs(ss**2, 1-cc**2)) == 0
    # 单点上 X_AA = X_BB（K11、K22 只经 K11K22 进入 Z）。
    assert sp.simplify(XAA-XBB) == 0 and sp.expand(Z.subs({k11: k11*zz, k22: k22/zz}, simultaneous=True)-Z) == 0
    vals = {}
    for (rv, sv) in [(1, sp.Rational(1, 3600)), (sp.Rational(1, 4), sp.Rational(1, 3600))]:
        av, h = sp.sqrt(rv), sp.sqrt(sv)
        for name, c_ in [("0", 1), ("pi/2", 0), ("pi", -1)]:
            v1 = sp.simplify(num.subs({a: av, hh: h, cc: c_, ss: sp.sqrt(1-c_**2)}))
            v2 = sp.simplify(XAB_closed.subs({R: rv, mu: h*c_*av**3, sg: sv}))
            assert v1.is_Rational and v1 == v2
            vals[rv, sv, name] = v1
    assert vals[1, sp.Rational(1, 3600), "0"] == sp.Rational(2747, 2058)
    assert vals[1, sp.Rational(1, 3600), "pi"] == sp.Rational(380123, 321602)
    assert vals[1, sp.Rational(1, 3600), "pi/2"] == sp.Rational(448253, 354482)
    assert vals[sp.Rational(1, 4), sp.Rational(1, 3600), "pi/2"] == sp.Rational(37114657, 31899904)
    assert vals[sp.Rational(1, 4), sp.Rational(1, 3600), "0"] == sp.Rational(43769977, 35617024)
    assert vals[sp.Rational(1, 4), sp.Rational(1, 3600), "pi"] == sp.Rational(1137131, 1051392)
    for rv in (1, sp.Rational(1, 4)):
        v0, vh, vp = (vals[rv, sp.Rational(1, 3600), k] for k in ("0", "pi/2", "pi"))
        assert v0 > vh > vp
    assert vals[1, sp.Rational(1, 3600), "0"]-vals[1, sp.Rational(1, 3600), "pi"] > sp.Rational(15, 100)
    # 正交支：a=0 时 mu=0，X 与 Phi 无关。
    assert sp.diff(XAB_closed.subs(back).subs(a, 0), cc) == 0
    # 首阶：X_AB = X_st(rho) + d_mu X|_0 * mu + O(sigma)，rho=1 处 d_mu X|_0 = 65/14。
    dmu = sp.diff(XAB_closed, mu).subs({mu: 0, sg: 0})
    assert sp.simplify(dmu.subs(R, 1)-sp.Rational(65, 14)) == 0
    xst = XAB_closed.subs({mu: 0, sg: 0})
    assert sp.simplify(xst.subs(R, 1)-sp.Rational(5, 4)) == 0
    lin = sp.simplify(dmu/xst - (120/(2+27*R+45*R**2+10*R**3)+240*(4+3*R)/p.subs(rho, R)-480/(1+18*R+45*R**2+20*R**3)))
    assert lin == 0
    return Z, vals


def main():
    t0 = time.time()
    op = charges()
    full_model_lattice()
    Z, vals = extended_site()
    print("耦合荷 (U(1)_Q, U(1)_L)：", op, "；eta W^3 的荷 (9,3)；唯一相位不变量 arg eta - 3 arg K21")
    print("完整模型：不变量 n arg(Fc)+m arg(eta_I)+... 须 n+3m=0，(3,-1) 有整数实现；对角汤川时无含 arg(Fc) 的不变量；去掉 Weinberg 结论不变")
    print("单点扩充探针 Z = 2211840[s^4 p(rho) + 120(eta K12^3 + c.c.)(4s+3K12K21) + 2880|eta|^2(s+3K12K21)]")
    for k, v in vals.items():
        print(f"  rho={k[0]}, sigma={k[1]}, Phi={k[2]}: X_AB = {v} = {float(v):.6f}")
    print("rho=1 处 dX_AB/dmu|_0 = 65/14；正交支与 Phi 无关")
    print(f"用时 {time.time()-t0:.1f} s")


if __name__ == "__main__":
    main()
