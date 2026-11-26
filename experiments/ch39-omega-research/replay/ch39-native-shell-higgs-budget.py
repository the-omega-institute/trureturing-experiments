"""原生十二个规范生成元和一个双重态，十九点单位边壳的一圈匹配预算。

设 R=theta(804)/theta(12)，r=R M_Z/mu。单个复标量双重态对
C=A1-3A2+2A3 的有限贡献为 F(lambda)/(5 pi)，其中
F(lambda)=(3/2) int_0^1 y^2 log[(lambda+aR^2)/(lambda+ar^2)] dy，
a=(1-y^2)/4。零模 F(0)=log(mu/M_Z)，标准模型 beta 已包含这一项。

对任意边相位 U，以 d 为最大度，有 |d I-L_U|=d I-L_0（逐元）。
幂级数遂给出 |exp(-t L_U)|<=exp(-t L_0)。由于上述对数比等于
int_0^infty exp(-t lambda) [exp(-tar^2)-exp(-taR^2)] dt/t，
F 是正的热核混合。整个壳替代原有轻标量时，Tr F(L_U)-F(0)
不超过平凡连接的十八个正模的贡献。有限 r>0 使零模减法有定义。

若保留外部轻 Higgs，只取 L_U 的较高十八模，需另作估计。删去极点
得到 Dirichlet 主子阵 D_U，交错性给出 lambda_(j+1)(L_U)>=lambda_j(D_U)。
同样的逐元反磁不等式适用于 D_U，故十八模曝光 <= Tr F(D_0)。
尖锐阈值下也成立：十八个高模均在参考与匹配尺度之间，且
log det D_U - log det D_0 = int [Tr exp(-t D_0)-Tr exp(-t D_U)] dt/t >= 0。
这些反磁比较针对原生质量矩阵，保持相同的非负裸质量也保留比较。
"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_mat

matter = run_path(str(Path(__file__).with_name("ch39-charge-matter-determinant-order.py")))
soft = matter["soft"]


def integer_spectrum(matrix):
    return sorted((z.real for z, multiplicity in matrix.charpoly().complex_roots()
                   for _ in range(multiplicity)), key=lambda v: v.mid())


def smooth(lambda_, momentum2):
    # (3/2) int_0^1 y^2 log(1+momentum2*(1-y^2)/(4 lambda_)) dy.
    k = momentum2/(4*lambda_+momentum2)
    v = k.sqrt()
    return (v.atanh()/v-1)/k-arb(1)/3


def budgets(levels, radius, reference_upper):
    assert all(reference_upper**2 < lam < radius**2 for lam in levels)
    raw = sum(smooth(lam, radius**2) for lam in levels)
    # log(1+x)<=x and (3/8) int_0^1 y^2(1-y^2)dy=1/20.
    reference_error = sum(reference_upper**2/(20*lam) for lam in levels)
    finite = raw-arb(0).union(reference_error)
    sharp = sum((radius/lam.sqrt()).log() for lam in levels)
    return finite/(5*arb.pi()), sharp/(5*arb.pi())


def verify():
    ctx.prec = 256
    h = acb.zeta_zero(11431).imag
    base, neutral, _ = soft["vector_pencil"](h)
    flat = base+neutral
    pole = soft["pole"]
    indices = list(range(19))
    integer = fmpz_mat([[int(flat[i,j].real.unique_fmpz()) for j in indices] for i in indices])
    keep = [i for i in indices if i != pole]
    dirichlet = fmpz_mat([[integer[i,j] for j in keep] for i in keep])
    assert all(integer[i,pole] == -1 for i in keep)  # star edges ensure D_U >= I
    assert max(integer[i,i] for i in indices) == 18  # Gershgorin: L_U <= 36 I
    flat_levels, dir_levels = integer_spectrum(integer), integer_spectrum(dirichlet)
    assert all(sum(integer[i,j] for j in indices) == 0 for i in indices)
    assert all(sum(dirichlet[i,j] for j in range(18)) == 1 for i in range(18))
    assert flat_levels[0].contains(0) and flat_levels[1] > 8 and dir_levels[0].contains(1)
    charged = base+matter["nonspokes"](3, h)
    levels = sorted((z.real for z in charged.eig()), key=lambda v: v.mid())
    assert len(levels) == 19 and arb("1.43399432") < levels[0] < arb("1.43399434")
    radius = soft["geo"]["theta"](804)/soft["geo"]["theta"](12)
    assert radius**2 > 36
    mz = arb("91.1876")
    reference_upper = radius*mz/arb("1e13")
    max_log = (arb("1e17")/mz).log()
    direction = (1, -3, 2)
    project = lambda v: sum(d*x for d, x in zip(direction, v))
    beta = (arb(41)/10, -arb(19)/6, arb(-7))
    beta_h = (arb(1)/10, arb(1)/6, arb(0))
    assert (project(beta)+arb(2)/5).contains(0)
    assert (project(beta)-project(beta_h)).contains(0) and project((0, 2, 3)) == 0
    em, weak = arb("127.930"), arb(".23129")
    initial = project((3*(1-weak)*em/5, weak*em, 1/arb(".1180")))
    center = initial+max_log/(5*arb.pi())
    assert center < arb("-10.60831")
    flat_budget = budgets(flat_levels[1:], radius, reference_upper)
    dir_budget = budgets(dir_levels, radius, reference_upper)
    upper18 = budgets(levels[1:], radius, reference_upper)
    entire = budgets(levels, radius, reference_upper)
    # A separately chosen common negative mass can tune the lowest eigenvalue
    # to zero. Recompute the actual q=3 gaps for that alternative explicitly.
    tuned = budgets([v-levels[0] for v in levels[1:]], radius, reference_upper)
    assert flat_budget[0] < arb(".22390")
    assert dir_budget[0] < arb(".27453") and dir_budget[1] < arb("1.01919")
    assert center+dir_budget[0] < arb("-10.33378") and center+dir_budget[1] < arb("-9.58912")
    assert upper18[0] < arb(".22600") and upper18[1] < arb(".93102")
    assert tuned[0] < arb(".24743") and tuned[1] < arb(".99464")
    print("Native unit graph, q=3, lowest squared mass / shell unit squared and mass / mu",
          levels[0].str(15), (levels[0].sqrt()/radius).str(15))
    print("SM central contrast at M_Z and maximum over 1e13--1e17 GeV", initial.str(15), center.str(15))
    print("Contributions to C, smooth finite subtraction / sharp activated logarithm")
    for label, values in (("flat eighteen positive modes", flat_budget),
                          ("any-phase eighteen unshifted high modes, upper bound", dir_budget),
                          ("q=3 eighteen unshifted high modes", upper18),
                          ("q=3 eighteen gaps above a separately tuned light mode", tuned),
                          ("q=3 entire nineteen-mode shell", entire)):
        print(label, *(v.str(15) for v in values))
    # An entire shell replacing the light Higgs has no separate light beta in C.
    print("q=3 full replacement, resulting C", *((initial+v).str(15) for v in entire))
    return flat_budget, dir_budget, upper18, tuned


if __name__ == "__main__":
    verify()
