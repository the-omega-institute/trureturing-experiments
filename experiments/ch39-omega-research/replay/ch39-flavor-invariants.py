"""双钟中性味模型内，由夸克质量平方和或 CKM 模方辨认六面相位。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx


flavor = run_path(str(Path(__file__).with_name("ch39-chiral-neutral-flavor-phase.py")))


def phase(c):
    return acb(0,c.acos()/3).exp()


def invariants():
    matrix = flavor["reference_matrix"]()
    mixing = flavor["ckm"]["standard_matrix"](*flavor["ckm"]["centers"])
    moduli = arb_mat([[abs(entry)**2 for entry in row] for row in mixing])
    nodes = [arb(i)/2 for i in range(-2,3)]
    data = []
    for c in nodes:
        y = matrix(phase(c))
        h = y*y.conjugate().transpose()
        characteristic = h.charpoly()
        masses = moduli.solve(arb_mat([[h[i,i].real] for i in range(3)]))
        discrepancy = sum((h*h)[i,i].real-masses[i,0]**2 for i in range(3))
        data.append([-characteristic[2].real,characteristic[1].real,-characteristic[0].real,discrepancy])
    polynomial = arb_mat([[c**i for i in range(5)] for c in nodes]).solve(arb_mat(data))
    # 各阶主子式由相应的复子式模方给出，次数上界见复算说明。
    return [arb_poly([polynomial[i,j] for i in range(degree+1)])
            for j,degree in enumerate((1,2,4,2))]


def reference_rotations():
    # 六面余弦 = Re(k U/W)。实部与虚部取自两组正交单位相位。
    real = flavor["singlet_at"]({724:acb(1),804:acb(1)})[1]
    imag = -flavor["singlet_at"]({724:acb(0,1),804:acb(1)})[1]
    k = acb(real,imag)
    c = flavor["singlet"](11661)[1]
    source = flavor["geom"]["phases"](acb.zeta_zero(8359).imag)
    relative = source[724]/source[804]
    return [(acb(c,sign*(abs(k)**2-c**2).sqrt())/(k*relative)).arg() for sign in (-1,1)]


def verify():
    ctx.prec=384
    c = flavor["singlet"](11661)[1]
    trace,second,third,discrepancy=invariants()
    other=-discrepancy[1]/discrepancy[2]-c
    assert arb("-147.405")<other<arb("-147.403")
    assert trace[1]<0 and discrepancy(c).contains(0)
    source_c=flavor["singlet"](8359)[1]
    trace_change=trace(source_c)/trace(c)-1
    assert arb("-7.580e-7")<trace_change<arb("-7.579e-7")
    angles=reference_rotations()
    assert arb("-.249218")<angles[0]<arb("-.249217")
    assert arb(".0500938")<angles[1]<arb(".0500939")
    print("CKM 必要二次式的另一根、质量平方和相对变化",other.str(14),trace_change.str(13))
    print("回到参考 CKM 的两个相对相位改变量",*(v.str(14) for v in angles))

    # 单个正欧氏质量参数 t 的有限费米子行列式。
    a,b,d=[p.derivative()(c) for p in (trace,second,third)]
    assert a<0<b and d>0
    t=(-b-(b*b-4*a*d).sqrt())/(2*a)
    determinant=t**3+t*t*trace+t*second+third
    assert determinant.derivative().derivative()(arb(0,1))<0
    assert determinant.derivative()(c).contains(0)
    momentum=t.sqrt()*arb("246.22")/arb(2).sqrt()
    assert arb(".64277912")<momentum<arb(".64277913")
    print("参考相位成为单模行列式唯一最大值时的欧氏动量 / GeV",momentum.str(14))


if __name__=="__main__":
    verify()
