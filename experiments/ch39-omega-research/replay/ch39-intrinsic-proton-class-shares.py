"""按壳面计数取权时，共同真空和编号零点的 Wilson 顶点实部下界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


# 四个驻点乘子和各竞争点的末项均除以 10^12，两个角坐标均除以 1000。
CERTIFICATES = (
    ((-1726364297, -579741649, -658139745, -335692209), (
        (-2478, 3126, 439981094), (-558, -634, 647336337),
        (-33, -2577, 693758762), (-25, 50, 6580802290),
        (-19, -2982, 8070309567), (-15, -150, 23116279019),
        (-14, -1176, 162002727), (-10, -150, 41411829632),
        (148, 2215, 2805093), (171, 613, 639974018), (205, -2092, 213045609))),
    ((-1120172922, -345936197, -539619121, -232951081), (
        (-338, -2811, 1628279824), (-146, -1855, 8398013),
        (-33, -2577, 181091132), (-30, -130, 1291682452),
        (-25, -125, 1855248139), (-19, -2982, 5614951045),
        (-15, -150, 42811384850), (-15, 30, 4285952734),
        (-14, -1176, 350732742), (147, -468, 267415801), (164, -44, 6801165481))),
    ((-3161945290, -934767286, -529164920, 347463614), (
        (-338, -2811, 2660407436), (-25, 50, 473132311),
        (-25, 55, 5246336357), (-19, -2982, 7070047926),
        (-15, -150, 95989196661), (-14, -1176, 338860086),
        (164, -44, 9633376418), (171, 613, 412997976))),
)
CLASSES = ((1, 0), (0, 1), (1, -1))
NUMBER = 2_306_888_100_085


def verify(root):
    def load(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    lib = load("mincut-global-winding")
    fixed = load("height-independent-flavor-vacuum")
    shell = load("three-ray-chiral-boundary")
    angles = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    ctx.prec = 224
    tri = lib["terms"](shell)
    sizes = [len(tri[k]) for k in CLASSES]
    assert sizes == [42, 42, 12] and len(tri[0, 0]) == 433
    assert all(m in (0, 1) for group in tri.values() for _, m in group)
    x = arb("-.01046")
    y = x-angles["reference_angles"](flavor)[0]
    height = acb.zeta_zero(NUMBER).imag
    phases = shell["geo"]["phases"](height)
    xn, yn = [phases[k].arg() for k in (724, 804)]
    shifts = (xn-x, yn-y, (xn-yn)-(x-y))
    deficit = arb(0).union(arb(23)/25*arb(2)**-20/14**2)
    duals = [([arb(i)/10**12 for i in mult],
              [(arb(a)/1000, arb(b)/1000, arb(lam)/10**12) for a, b, lam in cuts])
             for mult, cuts in CERTIFICATES]
    assert all(lam > 0 for _, cuts in duals for _, _, lam in cuts)
    # 单圈目标的剩余量按 (处方, 角位置) 收集，用精确正归一权线性组合。
    residuals = [[[], []] for _ in duals]
    inverse_amplitudes = []
    displacement_loss = []
    for q in range(1, 106):
        bare = lib["coeff"](tri, q, arb(0))
        softened = lib["coeff"](tri, q, deficit)
        largest = max(abs(z).lower() for z in softened)
        assert largest > 0
        inverse_amplitudes.append(1/largest)
        displacement_loss.append(q*q*sum(abs(c)*h*h for c, h in zip(softened, shifts))/192)
        grad = fixed["gradient_columns"](lib, tri, q, x, y)
        base = lib["energy"](q, bare, x, y)
        static = sum((1-m*deficit)*(z**q).real for z, m in tri[0, 0])
        constraints = []
        for grad_weights, cuts in duals:
            correction = sum(g*v for g, v in zip(grad_weights, grad))
            correction += sum(lam*(lib["energy"](q, bare, a, b)-base) for a, b, lam in cuts)
            constraints.append(correction)
        for j, (a, b) in enumerate(((x, y), (xn, yn))):
            z = [(c*acb(0, q*t).exp()).real for c, t in zip(softened, (a, b, a-b))]
            objectives = (sum(z)/96, sum(v/n for v, n in zip(z, sizes))/3,
                          (sum(z)+static)/529)
            for i, value in enumerate(objectives):
                residuals[i][j].append((value-constraints[i]).lower())
    fourier = 1/(192*sum(inverse_amplitudes))
    loss = max(v.upper() for v in displacement_loss)
    print("任意全局真空的傅里叶下界、指定位置的二阶损失上界", fourier.str(15), loss.str(15))
    assert fourier > arb(".00034702035") and loss < arb("3.628e-8")
    floors = [(min(ref), min(numbered)) for ref, numbered in residuals]
    for label, (ref, numbered) in zip(("96 面各自等权", "三个钟类各占三分之一", "全壳 529 面各自等权"), floors):
        print(label, "参考和编号位置的逐阶对偶下界", ref.str(15), numbered.str(15))
    # 96 面目标与负作用只差常数，所以目标在 A_* 的一阶导数消失。
    # 529 面目标仅增加一个角度常数，故二阶损失再乘 96/529。
    assert floors[0][0] > arb(".03223181") and floors[0][0]-loss > arb(".03223178")
    assert floors[1][0] > arb(".02562") and floors[1][1] > arb(".02557")
    assert floors[2][0] > arb(".07812720") and floors[2][0]-loss*96/529 > arb(".07812720")


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent)
