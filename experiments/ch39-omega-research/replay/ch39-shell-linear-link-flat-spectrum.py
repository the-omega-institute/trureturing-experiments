"""十九点壳的双基本边标量，扣除 Goldstone 后的平坦谱与可重整定向势。

每条边有一个复 5×5 标量。设 Φ_e=f_e U_e，f_e>0。一条边的单位取向
切空间为实 u(5)，在未破缺群下含色八重态、弱三重态、两个实单态和一个
复 (3,2)。规范轨道沿边为 ξ_i−Ad(U_e)ξ_j，ξ_i 属于 su(5)。在连通图上，
中性部分只有一份常量核，而带荷部分由非中心回路排除全部核。此处只在
物理实标量上计算指数，已被矢量吸收的 Goldstone 不再参加计数。

给定三角、四边迹系数时还能消除全部带荷角向零模。把辐边作为生成树，
以 D_j=Φ_oj/f_oj 记直达路径、L_ij=D_i Φ_ij/f_ij 记两步路径，反向边用
厄米共轭。在给定真空 D_j=1，L_ij=a_ij 1+b_ij P，P=diag(0,0,0,1,1)。
对每个非极点 j 选一个 b_rj≠0 的枢轴邻点 r，对其余邻点 i 加入

    Λ_ij⁴ ||L_ij − a_ij D_j − (b_ij/b_rj)(L_rj−a_rj D_j)||_F² .

正项均在目标真空为零，场次数至多四，展开的路径干涉项是三角、四边迹。
沿用原来的正径向势，零能方向只需在酉矩阵上研究。先用站点 U(5) 定树，
其 SU(5) 部分是规范自由度，其余为精确全局重相位。定树后这些等式要求
每个终点处 W_ij=a_ij 1+b_ij P_j。对任意 b≠0，酉性给出

    P_j† P_j + d P_j + conj(d) P_j† = 0,  d=conj(a)b/|b|² .

因 |a|=|a+b|=1，有 Re d=−1/2。两个不同的块相对相位 (a+b)/a 给出
不同的 Im d，相减迫使 P_j=P_j†，于是 P_j²=P_j。逆向边满足 W_ji=W_ij†，
故 b_ij≠0 的边使 P_i=P_j。在连通的非中心弦图上，所有 P_j 为同一投影。
同样相减线性化等式得 δP=δP† 和 PδP+δPP=δP，核仅含变换共同投影
的非对角方向，均为 SU(5) 共轭。平方和 Hessian 的核就是约束线性化的核。
因此目标给定秩二分支的全部物理角向零模只余18个站点相对迹重相位。
这里的系数逐边取自目标背景。驻点附近的谱隙与远处其他秩的真空需要区分。

单项场次数至多四时，仅用双基本边场和规范单态不能使用 SU(5) 的五指标
ε 构成非零站点净重相位的收缩，故该18维无质量空间在此类标量势下受保护。
"""

from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb


def shell_paths():
    geo = run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    pole = geo["sites"].index((17, 17, 17))
    clocks = geo["phases"](acb.zeta_zero(11431).imag)
    incoming = {j: {} for j in range(len(geo["sites"])) if j != pole}
    spokes = 0
    for i, j, distance, native in geo["edges"]:
        phase = native * clocks.get(distance, 1)
        if pole in (i, j):
            spokes += 1
            assert distance not in (200, 724, 804)
            assert geo["det"](geo["sites"][pole], geo["sites"][i], geo["sites"][j]) == 0
            assert phase.contains(1)
            continue
        first, last = phase**-2, phase**3
        incoming[j][i] = (first, last)
        incoming[i][j] = (first.conjugate(), last.conjugate())
        if not abs(phase**5 - 1) > arb("0.000001"):
            assert distance == 12
            assert geo["det"](geo["sites"][pole], geo["sites"][i], geo["sites"][j]) == 0
            assert phase.contains(1)
    assert (len(incoming), spokes, len(geo["edges"])) == (18, 18, 135)
    return incoming, len(geo["edges"])


def check_orientation(incoming):
    active = {}
    constraints = 0
    for j, paths in incoming.items():
        active[j] = {i for i, (a, z) in paths.items() if abs(z - a) > arb("0.000001")}
        ratios = {i: paths[i][1] / paths[i][0] for i in active[j]}
        r, t = next((r, t) for r, t in combinations(sorted(ratios), 2)
                    if abs(ratios[r] - ratios[t]) > arb("0.3987"))
        ar, zr = paths[r]
        at, zt = paths[t]
        dr = ar.conjugate() / (zr - ar).conjugate()
        dt = at.conjugate() / (zt - at).conjugate()
        assert (dr.real + arb(1) / 2).contains(0) and (dt.real + arb(1) / 2).contains(0)
        assert abs(dr.imag - dt.imag) > 0
        for i, (a, z) in paths.items():
            if i == r:
                continue
            factor = (z - a) / (zr - ar)
            direct = a - factor * ar
            assert (direct + factor * ar - a).contains(0)
            assert (direct + factor * zr - z).contains(0)
            constraints += 1
    noncentral = sum(map(len, active.values())) // 2
    assert (noncentral, min(map(len, active.values())), constraints) == (114, 9, 216)
    start = next(iter(active))
    seen, queue = {start}, [start]
    while queue:
        for i in active[queue.pop()]:
            if i not in seen:
                seen.add(i)
                queue.append(i)
    assert len(seen) == 18
    print("非中心弦、最小非中心度数、路径平方项", noncentral, 9, constraints)


def flat_spectrum(vertices, edges):
    neutral = edges - (vertices - 1)
    charged = edges - vertices
    singlets = edges + neutral
    goldstone = 24 * vertices - (8 + 3 + 1)
    physical = 8 * neutral + 3 * neutral + 12 * charged + singlets
    assert (neutral, charged, singlets, goldstone, physical) == (117, 116, 252, 444, 2931)
    assert physical + goldstone == 25 * edges
    # 一个复 (3,2) 的实化指数为 (5,3,2)。
    indices = (5 * charged, 3 * charged + 2 * neutral, 2 * charged + 3 * neutral)
    all_indices = tuple(t + 5 * edges for t in indices)  # 每边另有一个实径向伴随。
    assert indices == (580, 582, 583) and all_indices == (1255, 1257, 1258)
    print("中性伴随各份数、复 (3,2) 份数、实单态数、物理实角向数",
          neutral, charged, singlets, physical)
    print("平坦角向/全部物理边标量的实群指数", indices, all_indices)
    return indices


def check_running(indices):
    """有质量时，源逆耦合的欧氏一次减法给每个实标量贡献
    −T_i/(24π) ⟨log[(m²+x(1−x)Q²)/(m²+x(1−x)Q₀²)]⟩。
    Q≥Q₀ 时用 log(1+u)≤u，全部带荷标量的修正绝对值至多
    (Q²−Q₀²) Σ T_i/m² /(240π)。无质量项直接给 −T_i log(Q/Q₀)/(12π)。
    """
    # 在 M_Z 读入原中心值，仅跑标准模型与上述零质量角向，其他新增带荷态取重。
    weak, em, mass_z = arb(".23129"), arb("127.930"), arb("91.1876")
    initial = (3 * (1 - weak) * em / 5, weak * em, 1 / arb(".1180"))
    beta = tuple(base + Q(t, 6) for base, t in zip((Q(41, 10), Q(-19, 6), -7), indices))
    assert beta == (Q(3023, 30), Q(563, 6), Q(541, 6))
    poles = [mass_z * (2 * arb.pi() * a * b.denominator / b.numerator).exp()
             for a, b in zip(initial, beta)]
    assert arb("164.59079") < poles[2] < arb("164.59080") < poles[1] < poles[0]
    print("形式一圈逆耦合零点 Y/W/C，GeV", *(p.str(12) for p in poles))
    # 3(1−2x)² 为标量核的归一权，乘 x(1−x) 再积分控制重质量余项。
    weights = (3, -12, 12)
    assert sum(Q(c, i + 1) for i, c in enumerate(weights)) == 1
    assert sum(c * (Q(1, i + 2) - Q(1, i + 3)) for i, c in enumerate(weights)) == Q(1, 10)


if __name__ == "__main__":
    paths, edge_count = shell_paths()
    check_orientation(paths)
    check_running(flat_spectrum(len(paths) + 1, edge_count))
