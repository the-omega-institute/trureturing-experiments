"""可重整边标量的两个真空有相同树级质量和一圈规范阈值。

每条边用双基本矩阵 Φ_e，动能为 tr|DΦ_e|²。给定十九站点的壳背景
Φ_e=f_e U_e（f_e>0），另设一个具有正则动能的实规范单态 S。记

    F² = 10 ∑ f_e²,  n_e = √10 f_e/F,
    r_e = (tr Φ_e†Φ_e − 5 f_e²)/(√10 f_e),   R = ∑ n_e r_e,
    B_e = Φ_e†Φ_e − (tr Φ_e†Φ_e) 1₅/5.

取正数 b、λ_e、M、w，实数 c²+s²=1，并取规范不变势

    V = (b²/2) ∑ (r_e−n_e R)² + ∑ λ_e tr B_e²
        + (M²/2)(c R−s S)² + (w²/2)(s R+c S)².

每个单项的场次数至多四。给定壳背景和 S=0 使全部平方为零，是全局极小。
此例把边取向作为共同的平坦方向，选定真空只展示质量与顶点的可辨识性。
正则径向坐标 ξ_e 由 δΦ_e=ξ_e U_e/√10 给出，r_e=ξ_e+O(δΦ²)。
R 的线性部分 H=∑n_e ξ_e 是正则共同伸缩模。此势在正交径向和无迹
方向上给两模型相同的 Hessian。在 (H,S) 平面内，质量 w 的态为
h=s H+c S，另一态为 c H−s S，质量 M。

外部未破缺生成元逐边湮灭 U_e，故两个径向单态在外部规范背景中的
二次算子与背景无关。其余带荷标量及规范、规范固定与鬼场二次算子相同。
两模型的一圈二点读数由此相等，包含共同的额外边标量贡献。只计十维物质
和矢量塔的部分一圈读数与此处的完整一圈读数具有不同的场内容。

沿 H 有 Φ_e=f_e(1+H/F)U_e，矢量质量矩阵整体乘 (1+H/F)²。
故 h 对每个正则复矢量 X 的质量顶点为 G=2m_X² s/F，模指标对角。
外源的 X 反X 顶点也对角。对同一未破缺外源，A X h、A Z h 及 A X 反X h
顶点均为零。A X φ 对任一物理标量也为零，因为未破缺生成元保持规范
轨道切空间及其正交补。产生 X 反X h 的树图仅有两条矢量腿上的辐射。
下面用树质量和单位规范传播子算固定阶三体切割中的一个偏振振幅。
"""

from fractions import Fraction as Q


def check_mass_basis(cosine, sine, heavy, light):
    """质量矩阵写在 (H,S) 基底，两端检验归一的本征向量。"""
    assert cosine**2 + sine**2 == 1
    large = (cosine, -sine)
    small = (sine, cosine)
    matrix = [
        [heavy**2 * large[i] * large[j] + light**2 * small[i] * small[j]
         for j in range(2)]
        for i in range(2)
    ]
    for mass, direction in ((heavy, large), (light, small)):
        assert [sum(a * b for a, b in zip(row, direction)) for row in matrix] == [
            mass**2 * entry for entry in direction
        ]


def dot(a, b):
    return a[0] * b[0] - sum(x * y for x, y in zip(a[1:], b[1:]))


def difference(a, b):
    return tuple(x - y for x, y in zip(a, b))


def ym_vertex(momenta, polarizations):
    """三条腿全部取出射动量，省略共同的生成元矩阵元。"""
    return sum(
        dot(polarizations[i], polarizations[j])
        * dot(difference(momenta[i], momenta[j]), polarizations[k])
        for i, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1))
    )


def emission_diagrams(mass, scalar_mass, energy, momentum):
    """h 静止，带荷对沿 ±z 运动且均沿 x 极化，外源沿 z 极化。"""
    assert energy**2 - momentum**2 == mass**2
    charged = ((energy, 0, 0, momentum), (energy, 0, 0, -momentum))
    source = (-2 * energy - scalar_mass, 0, 0, 0)
    ex, ez = (0, 1, 0, 0), (0, 0, 0, 1)
    diagrams = []
    for leg in (0, 1):
        internal = (energy + scalar_mass, *charged[leg][1:])
        denominator = dot(internal, internal) - mass**2
        assert denominator == scalar_mass**2 + 2 * energy * scalar_mass > 0
        # 单位规范的极化投影，整体相位对两个图相同。
        projected = tuple(x - y * dot(internal, ex) / mass**2 for x, y in zip(ex, internal))
        vectors = [*charged, source]
        polarizations = [ex, ex, ez]
        vectors[leg], polarizations[leg] = internal, projected
        assert all(sum(row) == 0 for row in zip(*vectors))
        diagrams.append(ym_vertex(vectors, polarizations) / denominator)
    return diagrams


def verify():
    # 用最低带荷质量 m 为单位，振幅输出 F |𝓜|/(m |T|)。
    mass, light, heavy = Q(1), Q(1, 2), Q(7, 4)
    energy, momentum = Q(5, 4), Q(3, 4)
    threshold, sample = 2 * mass + light, 2 * energy + light
    assert (threshold, sample) == (Q(5, 2), 3)
    legs = emission_diagrams(mass, light, energy, momentum)
    assert legs == [Q(1), Q(1)]
    amplitudes = []
    for cosine, sine in ((Q(1), Q(0)), (Q(4, 5), Q(3, 5))):
        check_mass_basis(cosine, sine, heavy, light)
        amplitudes.append(abs(sum(legs) * 2 * mass * sine))
    assert amplitudes == [0, Q(12, 5)]

    # (3,2) 的六个复根，取 GUT 归一化超荷。
    indices = (6 * Q(25, 60), 3 * Q(1, 2), 2 * Q(1, 2))
    assert indices == (Q(5, 2), Q(3, 2), 1)
    print("三体阈值和样本质心能量，以 m 为单位", threshold, sample)
    print("两种势的偏振振幅 F |𝓜|/(m |T|)", *(str(value) for value in amplitudes))
    print("超荷、弱、色谱重的群指数", *(str(value) for value in indices))


if __name__ == "__main__":
    verify()
