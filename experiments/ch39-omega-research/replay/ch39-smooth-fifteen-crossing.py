"""单位边刚度下，把物质的阶跃激活换成两次欧氏动量之间的 Dirac 一圈减法。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_poly, ctx

partial = run_path(str(Path(__file__).with_name("ch39-partial-fifteen-crossing.py")))
matter, soft = partial["matter"], partial["soft"]


def upper_projection(blocks, radius):
    polynomials = {q:arb_poly([z.real for z in (-m).charpoly().coeffs()])
                   for q, m in blocks.items()}
    a, b, c = [polynomials[q] for q in (1, 4, 6)]
    # S+16C Dynkin coefficients by charge. A positive log derivative gives a deduction.
    numerator = -36*a.derivative()*b*c+arb("106.5")*b.derivative()*a*c-67*c.derivative()*a*b
    shifted = numerator(arb_poly([5, 1]))
    assert len(shifted.coeffs()) == 57 and all(v > 0 for v in shifted.coeffs())
    maximum = (11*b(0).log()-4*a(0).log()-7*c(0).log())/(5*arb.pi())
    assert maximum < arb(".17030")
    low_momentum = radius*arb("91.1876")/arb("1e13")
    for m in blocks.values():
        assert soft["lowest"](m, arb(0), arb(4)) > 1
    mean = 5+arb(270)/19
    cost = 7*19/(6*arb.pi())*(arb(4)/5*(1+radius**2/(4*mean)).log()-low_momentum**2/5)
    assert cost > 16*maximum
    return cost


def smooth_readout(block, mass2, radius, reference):
    # At this witness u*radius^2/lambda < 0.26, including u=1/4.
    matrix = block+mass2*soft["identity"]
    inverse = matrix.inv()
    assert radius**2/(4*mass2) < arb(".26")
    power = soft["identity"]
    exposure = arb(0)
    order = 28
    for k in range(1, order+1):
        power = power*inverse
        moment = 6*arb(factorial(k+1)**2)/factorial(2*k+3)
        exposure += ((-1)**(k+1)*moment*power.trace().real
                     *(radius**(2*k)-reference**(2*k))/(2*k))
    return exposure+arb(0, 19*arb(".26")**(order+1)/(order+1))


def verify():
    ctx.prec = 224
    height = acb.zeta_zero(11431).imag
    geo = soft["geo"]
    radius = geo["theta"](804)/geo["theta"](12)
    base, neutral, vector = soft["vector_pencil"](height)
    blocks = {q:base+matter["nonspokes"](q,height) for q in (1, 4, 6)}
    rbox = radius*arb(".999").union(arb("1.001"))
    cost = upper_projection(blocks, rbox)
    # Use certified enclosures of the exact step-fit parameters.
    mu = arb("1.03725489070050e14").union(arb("1.03725489070052e14"))
    mass2 = arb("70.7426557689436").union(arb("70.7426557689438"))
    reference = radius*arb("91.1876")/mu
    h = {q:smooth_readout(m, mass2, radius, reference) for q,m in blocks.items()}
    lift = acb_mat([[arb(1)/19]*19 for _ in range(19)])
    v0 = 18*radius.log()-(base+neutral+lift).det().real.log()/2
    v5 = 19*radius.log()-(base+vector).det().real.log()/2
    weak, em = arb(".23129"), arb("127.930")
    initial = [3*(1-weak)*em/5, weak*em, 1/arb(".1180")]
    beta = [arb(41)/10, -arb(19)/6, -7]
    indices = ((arb(1)/10,arb(8)/5,arb(9)/5), (arb(3)/2,0,2), (1,arb(5)/2,0))
    readings = [s-b*(mu/arb("91.1876")).log()/(2*arb.pi())
                +7*(j*v0+k*v5)/(4*arb.pi())
                -1770*2*sum(t*h[q] for t,q in zip(row,(1,4,6)))/(3*arb.pi())
                for s,b,j,k,row in zip(initial,beta,(0,2,3),(5,3,2),indices)]
    print("Small mass S cost lower and smooth inverse couplings at the step witness",
          cost.lower().str(9), *[v.str(13) for v in readings])
    assert all(v < 0 for v in readings)


if __name__ == "__main__":
    verify()
