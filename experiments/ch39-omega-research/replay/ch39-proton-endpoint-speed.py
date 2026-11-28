"""在第177枚正零点的面积权处认证保持一阶相消所需的最小端点总变差速率。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def area(x,y):
    return (x.conjugate()*y).imag


def verify():
    shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))
    ctx.prec = 192
    h = acb.zeta_zero(177).imag
    geo,soft = shell["geo"],shell["soft"]
    laplacian = shell["laplacian"](5,h)
    green = laplacian.inv()
    derivative = acb_mat(19,19)
    for i,j,d,_ in geo["edges"]:
        if d in (724,804):
            speed = 5*(geo["theta"](d)-geo["theta"](244))
            derivative[i,j] = acb(0,speed)*laplacian[i,j]
            derivative[j,i] = derivative[i,j].conjugate()
    weight_derivative = 2*h*(arb(23)/25)*arb(2)**-20/(h**2+soft["point"]["study"]["offset"])**2
    pole = soft["pole"]
    for j in soft["star"]["orbit"]:
        derivative[pole,pole] += weight_derivative
        derivative[j,j] += weight_derivative
        derivative[pole,j] -= weight_derivative
        derivative[j,pole] -= weight_derivative
    slope = -green*derivative*green

    sites = (0,2,3,5,6,9,11,13)
    support = (3,5,13)
    z3,z5,z13 = [green[6,j] for j in support]
    areas = [area(z5,z13),area(z13,z3),area(z3,z5)]
    assert all(value < 0 for value in areas)
    weights = [value/sum(areas) for value in areas]
    drift = sum(w*slope[6,j] for w,j in zip(weights,support,strict=True))
    first,second = z3-z5,z3-z13
    determinant = area(first,second)
    assert arb("-.00003234") < determinant < arb("-.00003233")
    r5 = area(-drift,second)/determinant
    r13 = area(first,-drift)/determinant
    assert r5 > 0 and r13 > 0
    optimum = r5+r13
    assert arb("6.39459637076796") < optimum < arb("6.39459637076798")

    # 该二阶实系统的逆给出两条活动转移在对偶界中取等号。
    dual = acb((second.imag-first.imag)/determinant,(first.real-second.real)/determinant)
    row = {i:sum(green[i,j]*w for j,w in zip(support,weights,strict=True)) for i in sites}
    jumps = {("a",6,i):row[i]-row[6] for i in sites if i!=6}
    jumps.update({("b",j,i):green[6,i]-green[6,j] for j in support for i in sites if i!=j})
    active = {("b",5,3),("b",13,3)}
    assert len(jumps)==28
    inactive = [(dual.conjugate()*z).real for key,z in jumps.items() if key not in active]
    assert max(inactive) < arb(".247")
    print("面积权处固定端点导数、5 到 3 和 13 到 3 的速率、最小总速率",
          drift.str(14),r5.str(15),r13.str(15),optimum.str(17))
    print("两活动方向面积、单侧方向总数、其余方向的对偶值上界",
          determinant.str(12),len(jumps),max(x.upper() for x in inactive).str(12))


if __name__ == "__main__":
    verify()
