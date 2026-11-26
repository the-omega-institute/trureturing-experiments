"""Bound the splitting of a real SU(4) adjoint supplying the colour threshold."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb

triplet = run_path(str(Path(__file__).with_name("ch39-published-color-triplet.py")))
data, comp, q, ratios = (triplet[k] for k in ("data","comp","q","ratios"))


def flow(three, eight):
    source = triplet["scalar_flow"]() if three else triplet["baseline"]
    if not eight:
        return source
    beta = source["BETA"].copy()
    matrix = [row.copy() for row in source["B"]]
    # Real colour adjoint: S_3=C_3=3, with half a complex scalar's weight.
    beta[2] += q("1/2")
    matrix[2][2] += 21
    return {**source, "BETA":beta, "B":matrix}


def corner_residual(onsets, upper, time):
    state = data["corner"](1,upper)
    h3,h8 = map(Fraction,onsets)
    breaks = sorted({Fraction(0),min(h3,time),min(h8,time),time})
    for a,b in zip(breaks,breaks[1:]):
        count = int(7*(b-a))+1
        step = q(b-a)/count
        running = flow(a >= h3, a >= h8)
        for _ in range(count):
            state = comp["advance"](state,step,upper,flow=running)
    target = comp["model"]["target"](11431,q(time))
    return [1/(r*u)-target for r,u in zip(ratios,state)]


def exclusion():
    strongest = flow(True,True)
    loop = comp["loop"]
    for i in range(3):
        slope = -2*(strongest["BETA"][i]+2*sum(max(x,0) for x in strongest["B"][i])/loop)/loop
        start = 1/(data["centers"][i]+data["widths"][i])**2
        assert min(start,start+33*slope) > arb(1)/2 and 1+slope > 2*comp["tolerance"]*ratios[i]
    boxes = [
        ((0,0),(33,25),"26.817553274557106"),
        ((0,25),("24.1",33),"26.817316742966472"),
        (("24.75",25),(33,33),"26.817173790274882"),
        ((0,25),("24.75",26),"26.81733976735785"),
        (("24.70",26),(33,33),"26.81718760855639"),
        (("24.1",26),("24.71","26.11"),"26.817277011003398"),
    ]
    for low,high,time in boxes:
        # Earlier thresholds lower the inverse-gauge vector field everywhere.
        lo = corner_residual(low,False,Fraction(time))
        hi = corner_residual(high,True,Fraction(time))
        assert max(lo) > comp["tolerance"] and min(hi) < -comp["tolerance"]
    assert (q("26.11")-q("24.70")).exp() > q("4.09")


def verify():
    exclusion()
    scale = q("26.81719035187015")
    h3, h8 = q("24.66780254876974"), q("26.21629520569882")
    values = data["centers"].copy()
    for i,sign in enumerate((1,-1,-1)):
        values[i] += sign*data["widths"][i]
    intervals = ((arb(0),h3,160,flow(False,False)), (h3,h8,20,flow(True,False)),
                 (h8,scale,12,flow(True,True)), (scale,arb(33),64,flow(True,True)))
    triplet["verify_ground"](values,intervals,2,scale,q("9.07e-8"),q(".136"))
    ratio = (h8-h3).exp()
    assert q("4.7043") < ratio < q("4.7044")
    print("required octet/triplet mass ratio above 4.09; witness",ratio.str(12))
    print("triplet and octet masses in GeV", (1000*h3.exp()).str(12),(1000*h8.exp()).str(12))


if __name__ == "__main__":
    verify()
