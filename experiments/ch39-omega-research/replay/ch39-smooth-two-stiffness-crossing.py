"""矢量与物质各用一套刚度时，两类平滑阈值的有限交点证书。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_mat

matter = run_path(str(Path(__file__).with_name("ch39-charge-matter-determinant-order.py")))
soft = matter["soft"]
indices = ((arb(1)/10, arb(8)/5, arb(9)/5), (arb(3)/2, 0, 2), (1, arb(5)/2, 0))
beta = (arb(41)/10, -arb(19)/6, -7)
count = 156


def masses(matrix):
    eigen = sorted((z.real for z in matrix.eig(multiple=True)), key=float)
    assert all(x > 0 for x in eigen)
    return eigen


def sheet(m2, q2, vector):
    k = q2/(4*m2+q2)
    a = k.sqrt().atanh()/k.sqrt()-1
    j0 = 2*a
    j2 = 2*(a/k-arb(1)/3)/3
    return (8*j0-3*j2)/14 if vector else 3*(j0-j2)/4


def setup():
    geo = soft["geo"]
    height = acb.zeta_zero(11431).imag
    star, neutral, heavy = soft["vector_pencil"](height)
    r = geo["theta"](804)/geo["theta"](12)
    eps = arb(1)/10**6
    scaled = star*10**6+neutral
    integer = fmpz_mat([[int(scaled[i,j].real.unique_fmpz()) for j in range(19)] for i in range(19)])
    neutral_masses = [z.real/10**6 for z, k in integer.charpoly().complex_roots()
                      if z != 0 for _ in range(k)]
    assert len(neutral_masses) == 18
    towers = [neutral_masses, masses(star+eps*heavy)]
    towers += [masses(star+20*matter["nonspokes"](q, height)) for q in (1, 4, 6)]
    assert min(v.lower() for seq in towers for v in seq) > arb("2e-6")
    return r, towers


def verify():
    ctx.prec = 224
    r, towers = setup()
    em, w = arb("127.930"), arb(".23129")
    sm = [3*(1-w)*em/5, w*em, 1/arb(".1180")]
    pi = arb.pi()
    def intercept(x, factor=1):
        q2 = (r*factor)**2
        def read(seq, mass, vector):
            return sum(sheet(lam+mass, q2, vector) for lam in seq)
        neutral, heavy = [read(seq, arb(0), True) for seq in towers[:2]]
        # Each exposure is strictly decreasing in its own positive bare mass.
        flavors = [read(seq, x.upper(), False).union(read(seq, x.lower(), False))
                   for seq in towers[2:]]
        return [bare+7*(n*neutral+b*heavy)/(4*pi)-2*count*sum(t*v for t,v in zip(row,flavors))/(3*pi)
                for bare,n,b,row in zip(sm,(0,2,3),(5,3,2),indices)]
    def differences(values):
        return values[0]-values[2], values[1]-values[2]
    def eliminate(values):
        first, second = differences(values)
        return (beta[1]-beta[2])*first-(beta[0]-beta[2])*second

    # Enclose the zero-reference root, then restore the actual M_Z subtraction
    # on all four faces. These two linear combinations have the same zeros.
    a, b = arb(0), arb(1)
    assert eliminate(intercept(a))*eliminate(intercept(b)) < 0
    positive_left = eliminate(intercept(a)) > 0
    for _ in range(55):
        mid = (a+b)/2
        sign = eliminate(intercept(mid))
        assert sign > 0 or sign < 0
        if (sign > 0) == positive_left:
            a = mid
        else:
            b = mid
    x = a.union(b)
    log_mu = 2*pi*differences(intercept(x))[0]/(beta[0]-beta[2])
    box = x+arb(0, arb("1e-9"))
    logs = log_mu+arb(0, arb("1e-5"))
    mz = arb("91.1876")
    # log(1+u*q0^2/lambda)/2 <= q0^2/(8*lambda) for either normalized kernel.
    q0 = r*(-logs).exp()
    per = q0**2/(8*arb("2e-6"))
    allowance = 19*per*(7*5/(4*pi)+2*count*arb("3.5")/(3*pi))
    assert allowance < arb("1e-17")
    error = arb(0, allowance.upper())
    def residual(mass, L, factor=1):
        return [t-b*(L+arb(factor).log())/(2*pi)+error
                for t,b in zip(intercept(mass,factor),beta)]
    e1 = eliminate(residual(box.lower(),logs.mid()))
    e2 = eliminate(residual(box.upper(),logs.mid()))
    # The light slope cancels exactly; only the bounded reference depends on L.
    assert e1*e2 < 0
    left = differences(residual(box,logs.lower()))[0]
    right = differences(residual(box,logs.upper()))[0]
    assert left > 0 and right < 0
    common = residual(box,logs)[2]
    assert arb("8.874e15") < mz*logs.exp() < arb("8.876e15")
    assert arb("6.6142") < common < arb("6.6144")
    assert arb(".0544042") < box < arb(".05440422")
    print("Smooth crossing: scale GeV, inverse, mass squared",
          (mz*logs.exp()).str(13), common.str(12), box.str(12))
    # Evaluate at a higher momentum holding the fitted masses fixed.
    later = residual(box,logs,arb("1.20"))
    print("Weak inverse at 1.20 times the crossing", later[1].str(11))
    assert later[1] < 0


if __name__ == "__main__":
    verify()
