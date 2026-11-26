"""Check the two bare ground sectors and a rotation-invariant contact that selects one."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb

joint = run_path(str(Path(__file__).with_name("ch39-joint-orientation-clock.py")))
fit = joint["fit"]
flavor = joint["contact"]["flavor"]


def verify():
    rows = joint["paired_scores"]()
    scores = {n:p for n,p,_,_ in rows}
    maximum = scores[11431]
    winners = []
    for order in (1,-1):
        for n in range(10815,12093):
            for first in (1,-1):
                for second in (1,-1):
                    if n==11431 and first==second==order:
                        winners.append((n,first,second,order))
                        continue
                    score = (scores[order*first*n]+scores[order*second*n])/2
                    assert maximum-score>arb(".00002295")
    assert len(winners)==2
    data = flavor["native"]
    p = data["p"]
    a,b = (data["geo"]["sites"][i] for i in data["portals"])
    def chord(x,y):
        return sum((u-v)**2 for u,v in zip(x,y,strict=True))
    assert (chord(p,a),chord(p,b),chord(a,b)) == (724,804,200)
    radius_sq = sum(x*x for x in p)
    swap_probability = F(chord(a,b),8*radius_sq)
    assert swap_probability == F(25,867)
    clock = joint["contact"]["clock"]
    bare = clock["free_ground"](maximum,2)
    _,d1 = clock["endpoint"](bare,1,2)
    _,d3 = clock["endpoint"](bare,3,2)
    endpoint_mass = 1/(4*(maximum*d1+(1-maximum)*d3))
    oriented_weight = next(w for n,p,w,obs in rows if n==11431)
    # The input rows already assign one half of the zero weight to each sign.
    probability = arb(swap_probability.numerator)/swap_probability.denominator
    bright_mass = endpoint_mass*oriented_weight/4*probability/625
    bare_gap = clock["free_ground"](maximum-arb(".00002295"),2)-bare
    assert 0<bare_gap<arb(".1") and bare+bare_gap < clock["free_ground"](arb(0),2)
    gap = bright_mass*bare_gap/(bare_gap+probability/2)
    assert gap > arb("4.14e-12")
    shifted = joint["contact"]["pulse_expectations"](
        joint["independent_copy_rows"](rows),joint["radial"]["local"]["graph_spectrum"](),
        2,strength=probability/2)[1]
    assert gap < shifted < bare_gap and arb("1.4118e-9") < shifted < arb("1.4119e-9")
    print("unscaled singlet contact: bright mass",bright_mass.str(9),"exact gap enclosed by",shifted.str(9))
    native = data["plane_response"](acb.zeta_zero(11431).imag)
    params = [arb("2.5605152"),arb("3.0832980")]
    def observables(reversed_fraction):
        effective = native+reversed_fraction*(native.conjugate()-native)
        w = fit["model"]["observables"]((arb(18)/25,*params),response=effective)
        return joint["contact"]["mixture"]["readouts"](joint["contact"]["mixture"]["direct_observables"](w))
    centers = flavor["shell"]["centers"]
    widths = [flavor["shell"]["sine_radius"]]*3+[flavor["shell"]["phase_radius"]]
    def fits(values):
        return all(abs(x-c)<w for x,c,w in zip(values,centers,widths,strict=True))
    assert fits(observables(arb(0)))
    # A fixed string of yes outcomes carries a coherent, rank-one boundary reference.
    # We use rational probabilities in the source and keep the full count before the first match.
    for ratio, expected in ((F(804,724),58),(F(3468-724,3468-804),203)):
        for count in range(expected+1):
            r = 1/(1+ratio**count)
            values = observables(arb(r.numerator)/r.denominator)
            if count==expected:
                assert fits(values)
            else:
                assert any(abs(x-c)>w for x,c,w in zip(values,centers,widths,strict=True))
        print("first passing comparison count", expected, "reverse population", float(r))
    print("bare doublet at 11431; source singlet probabilities 0 and 25/867")
    print("pure forward values",*(v.str(11) for v in observables(arb(0))))
    return maximum, endpoint_mass, probability


if __name__ == "__main__":
    verify()
