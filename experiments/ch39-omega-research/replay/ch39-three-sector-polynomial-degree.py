"""十三次的全类排除，含两个味部门及中性第11431枚零点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


PHASES = {
    -4: (
        (("0", "0"), ("1e-8", "0"), ("1e-8", "1e-8")),
        (("-1.1871596", ".1621108"), (".30485", ".5798161"),
         (".6343828", "-.9108542"), ("-1.1337387", ".2922111")),
    ),
    2: (
        (("0", "0"), ("-1e-8", "0"), ("-1e-8", "-1e-8")),
        ((".4221804", "-.7466266"), ("-.9437394", "-.035118")),
    ),
}


def verify():
    replay = Path(__file__).resolve().parent
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    selector = run_path(str(replay/"ch39-light-spectrum-zero-selector.py"))
    flavor = run_path(str(replay/"ch39-global-light-action-degree.py"))
    ctx.prec = 520
    degree = 13
    # 最接近原点的面很薄。由整数站位重算边相位，避免继承较低精度的包含半径。
    geo = shell["geo"]
    geo["edges"] = [(i,j,d,geo["unit"](i,j)*(-1 if d == 200 else 1))
                    for i,j,d,_ in geo["edges"]]
    theta = {d:(1-arb(d)/(2*geo["R"])).acos() for d in (244,724,804)}
    geo["phases"] = lambda h: {d:acb(0,(theta[d]-theta[244])*h).exp() for d in (724,804)}
    study = shell["soft"]["point"]["study"]
    study["offset"] = study["gauge"]["q"](study["gauge"]["OFFSET"])
    vectors = flavor["flavor_vectors"](shell,pump,selector,PHASES,degree)
    block = list(pump["rotation_bases"](shell))[1]
    adj = pump["adjoint"]

    def energy(n):
        matrix = shell["laplacian"](0,acb.zeta_zero(n).imag)
        return flavor["lowest"](adj(block)*matrix*block)

    lower,middle,upper = [energy(n) for n in (11430,11431,11432)]
    assert lower < middle < upper
    # 三角矩阵的相邻非对角项是关于能量的一阶差商，右上角为二阶差商。
    # 将 Chebyshev 递推作用于矩阵，避免直接相减两条高度极近的割线。
    slope = arb(2)/9
    nodes = arb_mat([[slope*lower-1,slope,0], [0,slope*middle-1,slope], [0,0,slope*upper-1]])
    previous = arb_mat([[1,0,0],[0,1,0],[0,0,1]])
    current = nodes
    following, combined = [], []
    for _ in range(degree):
        following.append(2*current[1,2])
        combined.append(2*current[0,2])
        previous,current = current,2*nodes*current-previous
    for v in (following,combined):
        length = sum(x*x for x in v).sqrt()
        assert length > 0
        vectors.append([x/length for x in v])
    radius = selector["interior_radius"](vectors)
    assert radius > arb("2.74e-14")
    print("十三维归一位移单形的原点面距下界",radius.str(15))


if __name__ == "__main__":
    verify()
