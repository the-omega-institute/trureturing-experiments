"""Keep both conjugate shell orientations in the static clock with equal prior weights."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

fit = run_path(str(Path(__file__).with_name("ch39-radial-gain-two-angle-fit.py")))
radial = fit["circle"]["two_ticks"]
contact = radial["contact"]


def paired_scores():
    early, late = radial["radial_reflection"]()
    rows = []
    for n, generator, weight, observables in contact["candidate_sources"](include_tangent=True):
        first, second = contact["pulses"]["port_outputs"](generator, 2, impulse=True)
        diagonal, other, real, imag = observables[5:9]
        for sign, immediate, delayed in zip((1,-1), first, second):
            probability = abs(early*immediate + late*delayed)**2
            data = [diagonal,other,real,sign*imag,arb(int(sign < 0)),observables[-1]]
            rows.append((sign*n,probability,weight,data))
    total = sum(weight for _,_,weight,_ in rows)
    return [(label,p,w/total,obs) for label,p,w,obs in rows]


def ground_readout(rows, modes):
    n,_,_,_,probability,_,means,residual = contact["pulse_expectations"](rows,modes,2)
    a,b,c,d,r,_ = means
    response = acb_mat([[a,acb(c,d)],[acb(c,-d),b]])
    assert n == 11431
    print("joint label", n, "probability", probability.str(11), "reversed mass", r.str(11))
    print("joint gauge bound", residual.upper().str(9))
    gain = arb(18)/25
    outer, params = fit["isolate"](("2.5605","3.0833"),response,gain)
    fit["prove_least"](response,gain,outer)
    # Also use the original two-angle fit to test whether it transfers without retuning.
    for pair in (params,[arb("2.5605152"),arb("3.0832980")]):
        w = fit["model"]["observables"]((gain,*pair),response=response)
        ud,ub,cd,cb = [w[i,j] for i,j in ((0,0),(0,2),(1,0),(1,2))]
        c13 = (1-abs(ub)**2).sqrt()
        phase = (-ud*ub.conjugate()*cd.conjugate()*cb-abs(ud*ub*cb)**2/(1-abs(ub)**2)).arg()
        values = [abs(w[0,1])/c13,abs(cb)/c13,abs(ub),phase]
        shell = fit["model"]["shell"]
        limits = [shell["sine_radius"]]*3+[shell["phase_radius"]]
        assert all(abs(v-c)<limit for v,c,limit in zip(values,shell["centers"],limits,strict=True))
        print("contact readouts",*(v.str(10) for v in values))
    return probability, r, residual


def direction_switch(rows):
    # Reverse the order of the source and detector with probability v, recording the choice.
    scores = {n:p for n,p,_,_ in rows}
    def lines(n):
        return scores[n], scores[-n]-scores[n]
    a,b = lines(11431)
    x,y = lines(11831)
    crossing = (a-x)/(y-b)
    assert arb(".24425162") < crossing < arb(".24425165")
    for n in scores:
        if n in (11431,11831):
            continue
        c,d = lines(n)
        assert a-c > 0 and a-c+(b-d)*crossing > 0
        assert x-c+(y-d)*arb(".25") > 0
    assert y>b
    print("first change of the maximum transmission", crossing.str(12), "11431 to 11831")


def thermal_bound(rows):
    peak = next(row for row in rows if row[0] == 11431)
    runner = max(p.upper() for n,p,_,_ in rows if n != 11431)
    clock = contact["clock"]
    bare = clock["free_ground"](peak[1],2)
    upper = bare+peak[2]/625
    gap = clock["free_ground"](runner,2)-upper
    assert gap > arb(".00004")
    assert clock["free_ground"](arb(0),2)-upper > gap and arb(".1") > gap+peak[2]/625
    excited_weight = 2**36*(-2**20*gap).exp()
    assert excited_weight < arb("1e-7")
    print("joint clock gap >", gap.lower().str(8), "Gibbs excited weight <", excited_weight.upper().str(8))


def independent_copy_rows(rows):
    # A separate fair orientation in each of two copies. The mixed signs have
    # identical scores, so their spectral weights add in the secular equation.
    paired = []
    for plus,minus in zip(rows[::2],rows[1::2],strict=True):
        n,forward,weight,pos = plus
        _,backward,_,neg = minus
        both = [(x+y)/2 for x,y in zip(pos,neg,strict=True)]
        paired.extend(((n,forward,weight/2,pos),(-n,backward,weight/2,neg),
                       (n+20000,(forward+backward)/2,weight,both)))
    return paired


def verify():
    rows = paired_scores()
    forward = max((row for row in rows if row[0]>0), key=lambda r:float(r[1].mid()))
    reverse = max((row for row in rows if row[0]<0), key=lambda r:float(r[1].mid()))
    assert forward[0]==11431 and reverse[0]==-10965
    assert all(reverse[1]>p for n,p,w,obs in rows if n<0 and n!=reverse[0])
    assert forward[1]-reverse[1]>arb(".000126")
    direction_switch(rows)
    modes = radial["local"]["graph_spectrum"]()
    p, r, error = ground_readout(rows, modes)
    assert p>arb(".999941") and r<arb("8.711e-6") and error<arb("1.805e-7")
    paired = independent_copy_rows(rows)
    p, r, error = ground_readout(paired, modes)
    assert p>arb(".999983") and r<arb("2.631e-6") and error<arb("1.216e-7")
    thermal_bound(paired)
    print("independent orientations:",p.str(12),r.str(12),error.upper().str(9))
    return rows


if __name__ == "__main__":
    verify()
