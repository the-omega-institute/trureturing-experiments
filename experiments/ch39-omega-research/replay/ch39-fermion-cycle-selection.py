"""四循环弱双态的一圈取向能差和四个质量阈值。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
from sympy import Matrix, diag, eye, factor, prod, symbols

winding = run_path(str(Path(__file__).with_name("ch39-first-arithmetic-plaquette-winding.py")))
angle = winding["angle"]


def determinant_check():
    # z 的形式共轭为 1/z。非等幅例与任意幅度的常数项分别检验。
    z, t = symbols("z t", nonzero=True)
    a = symbols("a:4")
    b = symbols("b:4")
    cycle = Matrix([[0, 0, 0, z], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]])
    matrix = diag(*a) + diag(*b)*cycle
    assert factor(matrix.det() - (prod(a)-prod(b)*z)) == 0
    numeric = matrix.subs(dict(zip((*a, *b), (2, 3, 5, 7, 1, 4, 6, 8))))
    gram = numeric.T.subs(z, 1/z) * numeric
    determinant = (t*eye(4)+gram).det(method="berkowitz")
    assert factor(determinant-determinant.subs(z, 1)+prod((2,3,5,7,1,4,6,8))*(z+1/z-2)) == 0


def masses_squared(phi):
    # M=1，r=1/2。每个数代表一份对角弱双态。
    return sorted([arb(5)/4 + (phi+k*arb.pi()/2).cos() for k in range(4)], key=float)


def quartic_energy(masses):
    # 等相位的项已扣去。两份弱分量，每份为四维 Dirac 场。
    return -2*sum(x**2*x.log() for x in masses)/(16*arb.pi()**2)


def shape(masses):
    return (masses[2]-masses[1])/(masses[3]-masses[0])


def compatible_shape(masses, relative_error):
    # 真实模型和待排除模型都可在相同测量值的误差窗中移动。
    ratio = ((1+relative_error)/(1-relative_error))**2
    compatible = [(x/ratio).union(x*ratio) for x in masses]
    return shape(compatible)


def verify():
    determinant_check()
    ctx.prec = 150
    # 每条循环边共用同一次 w。第一正零点已经处在所用单调弧内。
    assert 4*angle(1) < arb.pi()/6
    alpha = 2*arb.pi()/3
    direct = {}
    for sign in (-1, 1):
        for n in (11430, 11431, 11432):
            masses = masses_squared(angle(n) + sign*alpha)
            assert shape(masses).overlaps((arb.pi()/6 + sign*angle(n)).tan())
            direct[n,sign] = shape(masses)
            if n == 11431:
                error_box = compatible_shape(masses, arb("5e-10"))
        neighbours = sorted([direct[n,sign] for n in (11430,11432)], key=float)
        assert neighbours[0] < error_box < neighbours[1]
    assert direct[11431,-1] < 1/arb(3).sqrt() < direct[11431,1]

    # 其余三次循环跃迁为实，末次跃迁带 -conj(w)^p。
    power = 23226
    winding["verify"]()
    def seam(n, sign):
        return masses_squared((arb.pi()-power*angle(n))/4+sign*alpha)
    best, second = seam(11431, 1), seam(1144, 1)
    gap = quartic_energy(second)-quartic_energy(best)
    assert arb("6.4303e-13") < gap < arb("6.4305e-13")
    measurement = compatible_shape(best, arb("1e-6"))
    assert measurement.lower() > arb(".99997278") > arb(".99996948") > shape(second)
    assert quartic_energy(best) < quartic_energy(seam(11431, -1))
    # 将边界符号倒转，固定同一批相位时这两个候选的能量顺序反转。
    def unsigned(n):
        return masses_squared(-power*angle(n)/4 + alpha)
    assert quartic_energy(unsigned(11431)) > quartic_energy(unsigned(1144))
    print("四循环能差 / M^4", gap.str(15), "目标与第二名的质量读数",
          shape(best).str(15), shape(second).str(15),
          "两侧各容许 1 ppm 后的目标读数下界", measurement.lower().str(15))


if __name__ == "__main__":
    verify()
