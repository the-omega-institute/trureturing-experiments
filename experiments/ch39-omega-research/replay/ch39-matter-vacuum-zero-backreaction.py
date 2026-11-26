"""将全壳 Dirac 真空能反馈给每一枚实际正零点的味能量。"""
from collections import Counter
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

phase = run_path(str(Path(__file__).with_name("ch39-matter-number-three-spokes.py")))
soft = phase["soft"]
star = soft["star"]
pole, orbit = soft["pole"], star["orbit"]


def decrease(h):
    return arb(23)/25*arb(2)**-20/(h*h+soft["point"]["study"]["offset"])


def flavor_energy(h):
    # Replace the old single-site diagonal shift by the three actual spokes.
    mat = soft["flavor_block"](h,arb(1))
    mat[0,0] = 0
    loss = decrease(h)
    for j in orbit:
        mat[pole,pole] -= loss
        mat[j,j] -= loss
        mat[pole,j] += loss
        mat[j,pole] = mat[pole,j].conjugate()
    return soft["lowest"](mat,-arb(".001"),arb(0))


def charged_vacuum(h, cosine):
    # Hold the archived common matter angle fixed between the two shells.
    base = soft["vector_pencil"](h)[0]
    w = 1-decrease(h)
    angle = acb(cosine.mid(),(1-cosine.mid()**2).sqrt())
    shift = arb(0,(arb(3).sqrt()*cosine.rad()/(1-cosine**2).sqrt()).upper())
    trace = arb(0)
    for q in (-10,-5,5,10):
        mat = base+phase["crossing"]["matter"]["nonspokes"](q,h)
        for j in orbit:
            mat[pole,pole] += w-1
            mat[j,j] += w-1
            mat[pole,j] = -w*angle
            mat[j,pole] = mat[pole,j].conjugate()
        masses = [m+shift for m in phase["smooth"]["masses"](mat)]
        trace += (18 if abs(q)==5 else 3)*sum(m*m*m.log() for m in masses)
    return -trace/(16*arb.pi()**2)


def verify():
    _,_,cosine,radial_weights,_ = phase["verify"]()
    edges = [(i,j) for i,j,_,_ in soft["geo"]["edges"]]
    degrees = Counter(i for e in edges for i in e)
    h1 = acb.zeta_zero(1).imag
    delta = decrease(h1)
    dlow = {i:d-(3 if i==pole else int(i in orbit))*delta for i,d in degrees.items()}
    assert max(degrees.values()) == 18 and min(dlow.values())>9
    bound = 75*4*18*sum(1/(dlow[i]*dlow[j]).sqrt() for i,j in edges)/(3*8*arb.pi()**2)
    six = 3*sum(radial_weights)*(1-cosine)
    assert six < arb(".26")
    safe = arb("1.5e-17")/(bound+six)
    assert safe > arb("7e-20")
    print("Absolute magnetic vacuum bound and safe backreaction lower endpoint",
          bound.upper().str(9),safe.lower().str(10))

    h = acb.zeta_zero(11431).imag
    center_flavor = flavor_energy(h)
    center_charged = charged_vacuum(h,cosine)
    for index, floor, upper in ((2019542,".028","7e-11"),(93917305,".001","7.6e-13")):
        other = acb.zeta_zero(index).imag
        gap = flavor_energy(other)-center_flavor
        assert gap > 0
        charged = center_charged-charged_vacuum(other,cosine)
        # The neutral tower and the same-shell zero-connection reference only
        # change by the three-spoke stiffness. The trace derivative on [0,36]
        # is bounded by 300, including the continuous endpoint at zero.
        error = (33+75)*19*300*4*abs(decrease(h)-decrease(other))/(16*arb.pi()**2)
        advantage = charged+arb(0,error.upper())
        threshold = gap/advantage
        assert advantage > arb(floor) and threshold < arb(upper)
        print("Witness zero, flavor cost, loop advantage, trial crossover",
              index,gap.str(11),advantage.str(9),threshold.str(9))


if __name__ == "__main__":
    verify()
