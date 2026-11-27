"""三点有符号截止测度不能同时给出 K=1、W=c 与对共同截止标度的一阶驻定（|c| ≤ 1e-15）。

模型与第十一章及 ch39-cutoff-scale-stationarity.py 相同：电弱对称原点十六重标量的完整可见谱
（#11661、#11934、#12207），零裸项、无隐单态。截止测度 μ = Σ_j w_j δ_{x_j}，x_j ∈ (0,1]，Σ w_j = 1。

一阶驻定下（ch39-cutoff-scale-stationarity 说明的定理二）：K_μ=1 ⇔ Σ w L = −1，W_μ=c ⇔ Σ w P = 4c，DK_μ=0 ⇔ Σ w Φ₀=0，
DW_μ=0 ⇔ Σ w Ψ₀=0。故五个条件等价于四个齐次条件
    Σ w (1+L) = Σ w (P−4c) = Σ w Φ₀ = Σ w Ψ₀ = 0
外加归一。本证书证明：对 |c| ≤ 1e-15，任何三个互异截止上这四个齐次条件只有零解。

记 t = log x，A = 1+L > 0，s = −Ψ₀/A（严格增），f = Φ₀/A，g = (P−4c)/A。三点有非零解当且仅当
二阶差商 f[s₁,s₂,s₃] 与 g[s₁,s₂,s₃] 同时为零。差商等于 ∫ (·)_{ss} M dσ，M 为非负帐篷核。
以 ' 记 d/dt，N₁ = f''s' − f's''，N₂ = g''s' − g's''，则 f_ss = N₁/s'³，g_ss = N₂/s'³。
（i）μ₀ = 4e-12：在 t ≤ T* = −22.5 上 N₂ + μ₀N₁ < 0，故 t₃ ≤ T* 的三点 (g+μ₀f)[…] < 0。
（ii）μ₁ = 1e-7：在 I' = [−26.875, −25.3125] 之外 N₂ + μ₁N₁ < 0；帐篷核的凹性把 I' 上的正部
     与 [−25.3125, T*] 上的负部比较，得 t₃ ≥ T* 的三点 (g+μ₁f)[…] < 0。
两种情形下 f[…] = 0 都推出 g[…] < 0，故四个条件不能同时成立。
"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).resolve().parent
previous = run_path(str(ROOT / "ch39-cutoff-scale-stationarity.py"))

C_BOX = arb("1e-15").union(arb("-1e-15"))   # 证书对 |c| ≤ 1e-15 一并成立
MU0, MU1 = arb("4e-12"), arb("1e-7")
T_STAR, T_A, T_B, T_DEEP = arb("-22.5"), arb("-26.875"), arb("-25.3125"), arb(-40)
# 全部端点与步长为二进有理数，区间端点比较精确。
STEP, CUT = arb(1) / 32, arb(3000)


def decreasing(lo, hi, fn):
    """fn 在正半轴上递减时，[lo,hi] 上的像。"""
    return fn(hi).union(fn(lo))


def scaled_jets(rows, x, ulo_of):
    """在 x 的区间上包住缩放后的二阶射流。

    φ_m = Φ_m/x = Σαu^m e^{-u}，ψ_m = Ψ_m/x² = Σg u^m e^{-u}，ℓ = L/x，p = P/x²，q = Q/x²，
    Q = (1/2π)Σ g z² E₁(u)。递推 DL = Φ₁/2π，D²L = Φ₂/2π，DP = P+Q，DQ = Ψ₂/2π。
    返回 ŝ = s/x²（及其 D、D²），N̂₁ = N₁/x³，N̂₂ᴾ，N̂₂ᶜ，使 N₂ = x⁴N̂₂ᴾ + x³N̂₂ᶜ。
    """
    pi = arb.pi()
    ph, ps = [arb(0)] * 3, [arb(0)] * 3
    ell = p = q = tail = arb(0)
    for g, a, z in rows:
        if z.is_zero():
            ph[0] += a
            ps[0] += g
            continue
        ulo, uhi = ulo_of(z)
        if ulo > CUT:
            # u ≥ u_lo ≥ 2 时 u^m e^{-u}（m ≤ 2）、uE₁、uE₂、u²E₁ 都不超过 u_lo² e^{-u_lo}。
            tail += (abs(a) + abs(g)) * ulo * ulo * (-ulo).exp()
            continue
        u = ulo.union(uhi)
        e = decreasing(ulo, uhi, lambda v: (-v).exp())
        e1 = decreasing(ulo, uhi, lambda v: v.expint(1))
        e2 = decreasing(ulo, uhi, lambda v: v.expint(2))
        powers = [e,
                  decreasing(ulo, uhi, lambda v: v * (-v).exp()) if ulo > 1 else u * e,
                  decreasing(ulo, uhi, lambda v: v * v * (-v).exp()) if ulo > 2 else u * u * e]
        for m in range(3):
            ph[m] += a * powers[m]
            ps[m] += g * powers[m]
        ell += a * u * e1
        p += g * u * e2
        q += g * u * u * e1
    T = tail.union(-tail)
    ph = [v + T for v in ph]
    ps = [v + T for v in ps]
    ell, p, q = (ell + T) / (2 * pi), (p + T) / (2 * pi), (q + T) / (2 * pi)
    A0, A1, A2 = 1 + x * ell, x * ph[1] / (2 * pi), x * ph[2] / (2 * pi)
    assert A0 > 0

    def quotient(u0, u1, u2):
        v = u0 / A0
        dv = (u1 - v * A1) / A0
        return v, dv, (u2 - 2 * dv * A1 - v * A2) / A0

    s = quotient(-ps[0], -(2 * ps[0] + ps[1]), -(4 * ps[0] + 3 * ps[1] + ps[2]))
    f = quotient(ph[0], ph[0] + ph[1], ph[0] + ph[1] + ph[2])
    gP = quotient(p, p + q, p + q + ps[2] / (2 * pi))
    n1 = f[2] * s[1] - f[1] * s[2]
    n2p = gP[2] * s[1] - gP[1] * s[2]
    # g 中的 −4c/A：D(−4c/A) = 4cA₁/A²，D²(−4c/A) = 4c(A₂A − 2A₁²)/A³，二者均为 x 乘有界量。
    gc1 = 4 * C_BOX * ph[1] / (2 * pi * A0 * A0)
    gc2 = 4 * C_BOX * (ph[2] * A0 - 2 * x * ph[1] * ph[1] / (2 * pi)) / (2 * pi * A0 * A0 * A0)
    n2c = gc2 * s[1] - gc1 * s[2]
    return s, n1, n2p, n2c, (f, gP)


def box(rows, tlo, thi):
    xlo, xhi = tlo.exp(), thi.exp()

    def ulo_of(z):
        return arb((z / xhi).lower()), arb((z / xlo).upper())
    return (xlo.union(xhi),) + scaled_jets(rows, xlo.union(xhi), ulo_of)


def deep_tail(rows):
    """x ∈ [0, e^{−40}]：全部带质量项落入尾界，N̂₁ ≈ −124/6，s' > 0，两种组合都为负。"""
    xd = T_DEEP.exp()
    x = arb(0).union(xd)
    zmin = min((z.lower() for _, _, z in rows if not z.is_zero()), key=float)
    assert arb(zmin) / xd > CUT

    def ulo_of(z):
        return arb((z / xd).lower()), None
    s, n1, n2p, n2c, _ = scaled_jets(rows, x, ulo_of)
    # 零质量主项：ŝ = 62，ŝ' = 124，N̂₁ = (1/6)(124) − (1/6)(248) = −124/6。
    assert arb("123.9") < s[1] < arb("124.1") and n1 < arb("-20.6") and n1 > arb("-20.7")
    for mu in (MU0, MU1):
        assert x * n2p + n2c + mu * n1 < 0
    return n1


def signs(rows):
    """[−40,0] 上自适应二分：ŝ' > 0；t ≤ T* 上 h₀ < 0；I' 之外 h₁ < 0。"""
    count = [0, 0]

    def check(tlo, thi, depth):
        count[0] += 1
        x, s, n1, n2p, n2c, _ = box(rows, tlo, thi)
        ok = s[1] > 0
        if thi <= T_STAR:
            ok = ok and (x * n2p + n2c + MU0 * n1 < 0)
        if not (tlo >= T_A and thi <= T_B):
            ok = ok and (x * n2p + n2c + MU1 * n1 < 0)
        if ok:
            return
        assert depth < 8, ("区间未能定号", tlo.str(10), thi.str(10))
        count[1] += 1
        mid = (tlo + thi) / 2
        check(tlo, mid, depth + 1)
        check(mid, thi, depth + 1)

    t = T_DEEP
    while t < 0:
        check(t, t + STEP, 0)
        t = t + STEP
    return count


def tent_inequality(rows, dt=arb(1) / 256):
    """U = ∫_{I'} ψ₁⁺ dσ 的上和，与 Lw = ∫_{[t_B,T*]} ψ₁⁻ (σ*−σ)/(σ*−σ_A) dσ 的下和。

    ψ₁ dσ = (N₂+μ₁N₁)/s'² dt = ĥ₁/(x ŝ'²) dt，ĥ₁ = xN̂₂ᴾ + N̂₂ᶜ + μ₁N̂₁，σ = x²ŝ。
    """
    def sigma(t):
        x, s, *_ = box(rows, t, t)
        return x * x * s[0]
    s_star = arb(sigma(T_STAR).lower())
    s_a = arb(sigma(T_A).lower())
    denom = arb((s_star - s_a).upper())
    upper, lower = arb(0), arb(0)
    t = T_A
    while t < T_B:
        x, s, n1, n2p, n2c, _ = box(rows, t, t + dt)
        dens = (x * n2p + n2c + MU1 * n1) / (x * s[1] * s[1])
        if dens.upper() > 0:
            upper += arb(dens.upper()) * dt
        t = t + dt
    t = T_B
    while t < T_STAR:
        x, s, n1, n2p, n2c, _ = box(rows, t, t + dt)
        dens = -(x * n2p + n2c + MU1 * n1) / (x * s[1] * s[1]) * (s_star - x * x * s[0]) / denom
        # 被积函数非负（signs 已证 [t_B,T*] 上 h₁<0，且 σ ≤ σ*）；下界取 max(·,0)。
        if dens.lower() > 0:
            lower += arb(dens.lower()) * dt
        t = t + dt
    assert upper < lower
    return upper, lower


def point(rows, t):
    """点值 (s, f, g)，c 取 1e-122，用于独立的差商抽查。"""
    pi = arb.pi()
    x = t.exp()
    F = S = L = P = arb(0)
    for g, a, z in rows:
        if z.is_zero():
            F += a
            S += g
            continue
        u = z / x
        e = (-u).exp()
        F += a * e
        S += g * e
        L += a * z * u.expint(1)
        P += g * z * u.expint(2)
    A = 1 + L / (2 * pi)
    return -x * x * S / A, x * F / A, (x * P / (2 * pi) - 4 * arb("1e-122")) / A


def divided(v, s):
    return ((v[2] - v[1]) / (s[2] - s[1]) - (v[1] - v[0]) / (s[1] - s[0])) / (s[2] - s[0])


def spot_check(rows):
    """独立抽查（不是证明的一部分）：若干三点上直接由点值算差商，核对 H_μ < 0。"""
    triples = [(-30, -14, -12), (-27, -13.3, -12.6), (-26, -25.5, -10), (-24, -20, -16),
               (-40, -13, -1), (-25, -24.5, -24), (-20, -14, -13.9)]
    for tri in triples:
        pts = [point(rows, arb(t)) for t in tri]
        s = [p[0] for p in pts]
        F2 = divided([p[1] for p in pts], s)
        G2 = divided([p[2] for p in pts], s)
        mu = MU0 if tri[2] <= -22.5 else MU1
        assert G2 + mu * F2 < 0, tri
    return len(triples)


def jet_identity_check(rows):
    """有限差分核对 s、f、g 的递推射流（取谱的中点，256 位）。"""
    mid = [(g, arb(a.mid()), arb(z.mid())) for g, a, z in rows]
    h = arb("1e-15")
    for t0 in (arb(-25), arb(-14), arb(-3)):
        vals = []
        for t in (t0 - h, t0, t0 + h):
            x, s, _, _, _, (f, gP) = box(mid, t, t)
            vals.append([[x * x * v for v in s], [x * v for v in f], [x * x * v for v in gP]])
        for c in range(3):
            d1 = (vals[2][c][0] - vals[0][c][0]) / (2 * h)
            d2 = (vals[2][c][0] - 2 * vals[1][c][0] + vals[0][c][0]) / (h * h)
            assert abs(d1 / vals[1][c][1] - 1) < arb("1e-10")
            assert abs(d2 / vals[1][c][2] - 1) < arb("1e-10")


def chart_identity_check(rows):
    """与“#11661 的三点同侧排除”一节的记号对照（取谱的中点，点值核对，不是证明的一部分）。

    W(·,·,·) 记行 [F, DF, D²F] 的三阶 Wronski 行列式，Ψ = Ψ₀，R = W + Ψ/8π + c(Φ/2π − K)。
    书面证明：K−1 = Φ₀/2π − A，R = (P−4c)/4 + cA，故
      M_Φ = W(Ψ, K−1, Φ₀) = W(A, Ψ₀, Φ₀) = −A³N₁，
      M_R = W(Ψ, K−1, R) = −(A³/4)N₂ + (1/8π)[W(Ψ₀, Φ₀, P−4c) + 4c W(Ψ₀, Φ₀, A)]。
    返回额外项与主项之比在所查各点上的最大值。
    """
    pi, c = arb.pi(), arb("1e-122")
    mid = [(g, arb(a.mid()), arb(z.mid())) for g, a, z in rows]

    def wr(a, b, d):
        return (a[0] * (b[1] * d[2] - b[2] * d[1]) - a[1] * (b[0] * d[2] - b[2] * d[0])
                + a[2] * (b[0] * d[1] - b[1] * d[0]))
    worst = arb(0)
    for t in (arb(-26), arb("-21.5"), arb("-14.5"), arb(-10), arb(-3)):
        x = t.exp()
        F = [arb(0)] * 3
        S = [arb(0)] * 3
        K = W = L = P = Q = arb(0)
        for g, a, z in mid:
            if z.is_zero():
                F[0] += a
                S[0] += g
                K += a
                W += arb(g) / 2
                continue
            u = z / x
            e = (-u).exp()
            for m in range(3):
                F[m] += a * u**m * e
                S[m] += g * u**m * e
            K += a * u.expint(2)
            W += g * u.expint(3)
            L += a * z * u.expint(1)
            P += g * z * u.expint(2)
            Q += g * z * z * u.expint(1)
        Phi = [x * v for v in F]
        Psi = [x * x * v for v in S]
        K, W = x * K / (2 * pi), -x * x * W / (4 * pi)
        L, P, Q = L / (2 * pi), x * P / (2 * pi), Q / (2 * pi)
        phi0 = (Phi[0], Phi[0] + Phi[1], Phi[0] + Phi[1] + Phi[2])
        psi0 = (Psi[0], 2 * Psi[0] + Psi[1], 4 * Psi[0] + 3 * Psi[1] + Psi[2])
        A = (1 + L, Phi[1] / (2 * pi), Phi[2] / (2 * pi))
        Pc = (P - 4 * c, P + Q, P + Q + Psi[2] / (2 * pi))
        # K、W 的射流直接由定义给出：DK = Φ₀/2π，DW = −Ψ₀/4π。
        km1 = (K - 1, phi0[0] / (2 * pi), phi0[1] / (2 * pi))
        Wj = (W, -psi0[0] / (4 * pi), -psi0[1] / (4 * pi))
        R = tuple(Wj[i] + psi0[i] / (8 * pi) + c * (phi0[i] / (2 * pi) - (K if i == 0 else km1[i]))
                  for i in range(3))
        assert abs(R[0] - (P / 4 + c * L)) < arb("1e-60") * (abs(P) + 1)
        M_phi = wr(psi0, km1, phi0)
        M_R = wr(psi0, km1, R)
        xb, s, n1, n2p, n2c, _ = box(mid, t, t)
        h1 = xb * n2p + n2c + MU1 * n1
        # 该节的三个反例点：本文组合在 −26（位于 I′ 内）为正，在 −21.5、−14.5 为负。
        if t == arb(-26):
            assert h1 > 0
        elif t in (arb("-21.5"), arb("-14.5")):
            assert h1 < 0
        N1 = x**3 * n1
        assert abs(M_phi - wr(A, psi0, phi0)) < arb("1e-50") * abs(M_phi)
        assert abs(M_phi + A[0]**3 * N1) < arb("1e-30") * abs(M_phi)
        main = -wr(psi0, A, Pc) / 4
        extra = (wr(psi0, phi0, Pc) + 4 * c * wr(psi0, phi0, A)) / (8 * pi)
        assert abs(M_R - main - extra) < arb("1e-50") * (abs(main) + abs(extra))
        # c = 0 部分与证书中的 N̂₂ᴾ 一致：W(Ψ₀, A, P) = A³ x⁴ N̂₂ᴾ。
        Pz = (P, P + Q, P + Q + Psi[2] / (2 * pi))
        assert abs(wr(psi0, A, Pz) - A[0]**3 * x**4 * n2p) < arb("1e-30") * abs(wr(psi0, A, Pz))
        worst = worst.max(abs(extra / main))
    return worst


EXTRA_BOUND = {11661: "1.2e-5", 11934: "1.2e-5", 12207: "1.2e-5"}
EXPECT = {11661: ("152.7", "152.8", "160.8", "160.9"), 11934: ("153.5", "153.6", "161.3", "161.4"),
          12207: ("152.5", "152.6", "160.5", "160.6")}


def verify():
    data = previous["spectra"]()
    ctx.prec = 128
    for n, rows in data.items():
        assert sum(g for g, _, _ in rows) == -310
        massless = [(g, a) for g, a, z in rows if z.is_zero()]
        assert len(massless) == 1 and massless[0][0] == -62 and (massless[0][1] - arb(1) / 6).contains(0)
        assert all(z > arb("1.5e-11") for _, _, z in rows if not z.is_zero())
        n1 = deep_tail(rows)
        boxes, splits = signs(rows)
        U, Lw = tent_inequality(rows)
        lo_u, hi_u, lo_l, hi_l = (arb(v) for v in EXPECT[n])
        assert lo_u < U < hi_u and lo_l < Lw < hi_l and Lw > arb("1.04") * U
        ctx.prec = 256
        jet_identity_check(rows)
        k = spot_check(rows)
        extra = chart_identity_check(rows)
        assert extra < arb(EXTRA_BOUND[n])
        ctx.prec = 128
        print(n, "尾段 N̂₁", n1.str(5), "；[−40,0] 定号区间", boxes, "（二分", splits, "次）")
        print("   帐篷不等式：I' 正部上和 U", U.str(6), "< 加权负部下和 Lw", Lw.str(6))
        print("   M_R 额外项与主项之比（中点谱，五点）的最大值", extra.str(3))
        print("   结论：|c| ≤ 1e-15 时三个互异截止上四个齐次条件只有零解；抽查三点", k, "组")


if __name__ == "__main__":
    verify()
