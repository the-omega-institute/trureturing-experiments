"""十九点壳上同一 U(1) 链的电荷幂、轻零模与全部高模。

取正边权，整数荷 q 的平方质量算子为 m0² I + s² L_q，s>0、m0²≥0，
其磁拉普拉斯的能量是 Σ w_ij |f_i−U_ij^q f_j|²。存在质量为零的荷一
物质，便有 m0=0，且沿每条边平行。连通性迫使所有 |f_i| 相等且非零，
任意闭路的输运乘积都为1。反之，定树平凡化即得到一维复零空间。
同一站点规选取 q 次幂后，L_q 与 L_0 幺正等价。任意一组须保留零模的
非零整数荷，等价条件是每条闭路的输运的阶数整除这些荷绝对值的最大公因数。
三代 Weyl 规范量子数含 q=6Y=1。这里每个内部标量拉普拉斯只给一份空间
零模；L_−q 与 L_q 共轭等谱，手征标签与家族多重性在此处是外部输入。

稳定性。令 ε=λ_min(L_1)，f 为归一化基态，r_i=|f_i|。取 γ 为 L_0 的
谱隙下界。由 |r_i−r_j|≤|f_i−U_ij f_j| 与 Poincaré 不等式，
有 Σ(r_i−mean(r))²≤ε/γ。19 个分量之和的正交投影给
min r_i≥[√(1−ε/γ)−√(18 ε/γ)]/√19=:m。m>0 时使基态各分量为正，
边能量分成 (r_i−r_j)²+r_i r_j |1−U'_ij|²。因此 Σw|1−U'|²≤ε/m²，
任意闭路 C 满足 |1−hol(C)|≤√(Σ_C 1/w) √ε/m（简单闭路）。
由 |1−z^q|≤|q| |1−z| 和矩阵 Frobenius 范数，排序后的任意指标 k 满足
|λ_k(L_q)−λ_k(L_0)|≤|q|√(2 max(w) ε)/m。同时，以 r_i (f_i/r_i)^q
作试探向量（r_i=0 时取零），逐边能量至多为 q² 倍原能量，因此对任意
ε 和 |q|≥1，λ_min(L_q)≤q² ε。维数与边权变化时相同推导须替换 n、γ。
"""

from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx, fmpz_mat


def integer_spectrum(matrix):
    integer = fmpz_mat([[int(matrix[i, j].real.unique_fmpz()) for j in range(19)]
                        for i in range(19)])
    return sorted([z.real for z, count in integer.charpoly().complex_roots()
                   for _ in range(count)], key=float)


def laplacian(edges, transport):
    matrix = acb_mat(19, 19)
    for i, j, link in edges:
        matrix[i, i] += 1
        matrix[j, j] += 1
        matrix[i, j] = -transport(i, j, link)
        matrix[j, i] = matrix[i, j].conjugate()
    return matrix


def check_shell():
    base = run_path(str(Path(__file__).with_name("ch39-native-weyl-graph-diameter.py")))
    soft = base["soft"]
    geo = soft["geo"]
    height = acb.zeta_zero(11431).imag
    phases = geo["phases"](height)
    pole = soft["pole"]
    edges = [(i, j, acb(1) if pole in (i, j) else z * phases.get(d, 1))
             for i, j, d, z in geo["edges"]]
    assert len(edges) == 135
    flat = integer_spectrum(laplacian(edges, lambda i, j, z: 1))
    assert flat[0].is_zero()
    decrement = arb(23) / 25 * arb(2)**-20 / (height**2 + soft["point"]["study"]["offset"])
    assert 0 < decrement < arb("1e-12")
    # 对每种相位均有 (1−decrement)L_unit ≤ L_soft ≤ L_unit。
    gap = (1 - decrement) * flat[1]
    assert gap > arb("8.24686")
    actual = laplacian(edges, lambda i, j, z: z)
    ground = min(z.real for z in actual.eig(multiple=True))
    assert arb("2.75157169079") < (1 - decrement) * ground
    assert ground < arb("2.75157169080")
    print("平坦谱隙下界、实际11431荷一底谱（微调前）", gap.lower().str(13), ground.str(15))

    # 弦 (0,2) 的符号可由穿过极点的三角闭路测出。偶荷全平坦，奇荷全同谱。
    marked = (0, 2)
    pairs = {(i, j) for i, j, _ in edges}
    assert all(tuple(sorted(edge)) in pairs for edge in (marked, (pole, 0), (pole, 2)))
    odd = integer_spectrum(laplacian(edges, lambda i, j, z: -1 if (i, j) == marked else 1))
    assert arb(".17631277866") < odd[0] < arb(".17631277868")
    print("单弦负号时，偶荷与奇荷的最低本征值", 0, odd[0].str(13))


def stability_and_matching():
    gap, epsilon = arb("8.24686"), arb("1e-6")
    variance = epsilon / gap
    radius = ((1 - variance).sqrt() - (18 * variance).sqrt()) / arb(19).sqrt()
    shift = (2 * epsilon).sqrt() / radius  # max(w) ≤ 1
    assert radius > arb(".22907") and shift < arb(".006174")
    # 三角上由逐边同一规选给 gauge-invariant 通量的弦长上界。
    triangle = (arb(3) / (1 - arb("1e-12")) * epsilon).sqrt() / radius
    assert triangle < arb(".007562")

    # 对 j≥2，所有平坦本征值≥gap，质量比与有共同非负裸质量时的比
    # 都介于 sqrt(1−q shift/gap) 和其倒数之间。单高模有限减法响应
    # 对 log m 为 1-Lipschitz，任意公共尺度均可直接加和这 18 个误差。
    def log_error(q):
        return -(1 - q * shift / gap).log() / 2

    pi = arb.pi()
    qabs = (1, 4, 2, 3, 6)
    weyl = (-172, 128, 84, -84, 44)
    assert sum(weyl) == 0
    correction = 18 * log_error(3) / (5 * pi) + 3 * 18 / (2 * pi) * sum(
        abs(w) * log_error(q) / 105 for w, q in zip(weyl, qabs))
    # 共同平坦谱时 Weyl 投影抵消。标量激活核不超过矢量核，使两者之和≤0。
    em, weak = arb("127.930"), arb(".23129")
    hyper, left, color = 3 * (1 - weak) * em / 5, weak * em, 1 / arb(".1180")
    log_min = (arb("1e13") / arb("91.1876")).log()
    baseline = hyper - 3 * left + 2 * color + (5 * hyper - 3 * left - 2 * color) / 105
    baseline -= log_min / (105 * pi)
    bound = baseline + correction
    assert correction < arb(".043256") and bound < arb("-11.0438")
    assert -bound * 35 / 106 > arb("3.6465")
    print("近零模的各荷全谱误差系数、三角通量弦长上界", shift.upper().str(12),
          triangle.upper().str(12))
    print("18高模投影的新增误差、P*上界、三支极差下界", correction.upper().str(12),
          bound.upper().str(12), (-bound * 35 / 106).lower().str(12))


if __name__ == "__main__":
    ctx.prec = 224
    check_shell()
    stability_and_matching()
