"""将三低模以外的十六层 Dirac 物质纳入共享边的一圈阈值。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    directory = Path(__file__).parent
    shell = run_path(str(directory/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 256
    height = acb.zeta_zero(11431).imag
    rows = [[2*x for x in row] for row in shell["model"]["weyl"]]
    charges = shell["model"]["charges"]
    corners = []
    for em in (arb("127.90"),arb("127.96")):
        for s in (arb(".23126"),arb(".23132")):
            a1, a2, a3 = 3*(1-s)*em/5, s*em, 1/arb(".076")
            b1,b2,b3 = arb(41)/10,-arb(19)/6,-7
            n = ((a2*b3-a3*b2)/(a1*b2-a2*b1), (a3*b1-a1*b3)/(a1*b2-a2*b1), arb(1))
            assert n[0] > arb(".536") and n[1] > -arb("1.516")
            coefficients = [sum(a*b for a,b in zip(n,row)) for row in rows]
            negative = -sum(c for c in coefficients if c < 0)
            assert all(c > 0 or c < 0 for c in coefficients)
            assert sum(coefficients) > 0 and negative < arb("2.423")
            corners.append(coefficients)
    mean = (arb(270)-6*arb(0).union(arb("1e-8")))/19
    variance = (arb(4290)/19-mean**2).upper()
    spread = arb(0)
    for j in range(4,20):
        low = mean.lower()-(variance*(19-j)/j).sqrt()
        high = min(arb(36), mean.upper()+(variance*(j-1)/(20-j)).sqrt())
        assert low > 0
        spread += (high/low).log()/2
    cost = arb("2.423")*spread/(2*arb.pi())
    bound = 1/(1/arb(".076")-cost)
    assert spread < arb("6.185") and cost < arb("2.386") and bound < arb(".09284")
    print("任意相位、十六层对数谱差上界、阈值预算、一圈强耦合上界",
          spread.str(11),cost.str(11),bound.str(11))
    # Below all knots the primitive is zero, above the last it increases.
    knots = []
    for i,q in enumerate(charges):
        eigen = sorted((z.real for z in shell["laplacian"](q,height).eig()),key=float)
        knots.extend((x.log()/2,i) for x in eigen[3:])
    knots.sort(key=lambda item:float(item[0]))
    assert all(x < y for (x,_),(y,_) in zip(knots,knots[1:]))
    losses = []
    for coefficients in corners:
        primitive = slope = peak = loss = arb(0)
        last = knots[0][0]
        for point,i in knots:
            primitive += slope*(point-last)
            peak = peak.max(primitive)
            loss = loss.max(peak-primitive)
            slope += coefficients[i]/(2*arb.pi())
            last = point
        assert slope > 0
        losses.append(loss)
    worst = max(x.upper() for x in losses)
    fixed_bound = 1/(1/arb(".076")-worst)
    assert worst < arb(".0463") and fixed_bound < arb(".07627")
    print("11431 的全标尺最大投影回落、带十六层后的强耦合上界",
          worst.str(12),fixed_bound.str(12))


if __name__ == "__main__":
    verify()
