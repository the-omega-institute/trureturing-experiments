"""Find the least allowed maximum in-box word length for independent port mixtures."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
floating = run_path(str(root/"ch39-floating-clock-ground.py"))
gauge = floating["gauge"]
shortest, events = gauge["shortest"], gauge["events"]
model, words, q = shortest["model"], shortest["words"], shortest["q"]
counts = (97,279,136)
coordinates = ((-235723,110294,150368,249494),(-526967,246408,335994,558244),
               (-527263,246666,336164,557748))


def exclude_ports(ranges):
    time = q("26.81676714125598")
    _,_,logs = events["mixing"]["volume"]["input_bounds"]()
    inverse = [q(lo).union(q(hi)) for lo,hi in ranges]
    end,_,weak_minus_color = model["evolve"](inverse,logs,time,8192)
    weak = end[1]/shortest["ratios"][1]-model["target"](11431,time)
    speeds = floating["flow_bounds"]()
    anchor = q("429/16")
    target = model["target"](11431,anchor+q("1/8")).union(model["target"](11431,anchor-q("1/8")))
    beta = 1+(speeds[2]-speeds[1])/(2*target+speeds[1])
    residual = (1-beta)*weak-weak_minus_color
    assert residual.lower() > q(".00013") > (1+beta)*(arb(2)**-21+arb("1.3e-16"))
    return residual.lower()


def mixed_ports():
    mass = [0]*4
    for num,coords,depth in zip(counts,coordinates,(16,17,17),strict=True):
        assert sum(v*v for v in coords) == 5**depth
        word = words["factor"](coords,depth)
        quat = (1,0,0,0)
        for j in word:
            quat = words["multiply"](quat,words["GENERATORS"][j])
        assert quat==coords and len(word)==depth
        assert all(lo < abs(Fraction(v,coords[0])) < hi
                   for v,(lo,hi) in zip(coords[1:],words["BANDS"],strict=True))
        for j,v in enumerate(coords):
            mass[j] += num*5**(17-depth)*v*v
    assert sum(mass)==512*5**17
    return mass


def common_state(ports):
    for copies in range(1,33):
        assert not events["possible_light_events"](events["orbit_weights_from_ports"](copies,ports))
    ys = table_readout(ports)
    assert all(events["low"][n]<r<events["high"][n] for n,r in zip(events["names"],ys,strict=True))
    return [arb(v)/ports[0] for v in ports[1:]]+[q(r) for r in ys]


def table_readout(ports):
    weights = events["orbit_weights_from_ports"](33,ports)
    masses,used = events["event_masses"](weights)
    extra = (18,6,3)
    table = (root/"ch39-common-yukawa-events.tsv").read_text().splitlines()
    assert extra in used and any(line.split() == ["18","6","3","ref"] for line in table)
    masses["ref"] -= weights[extra]
    return [Fraction(masses[n],masses["ref"]) for n in events["names"]]


def residuals(state,time):
    for _ in range(128):
        state = shortest["reference"]["advance"](state,time/128,14)
    target = model["target"](11431,time)
    mix = [a+b for a,b in events["mixing"]["displacement"](27)]
    return [abs(1/(r*u)-target)+d for r,u,d in zip(shortest["ratios"],state[:3],mix,strict=True)]


def verify():
    assert all(not words["hits"](d) for d in range(16))
    coords = [v for v,w in words["hits"](16)]
    assert len(coords)==4
    ranges = []
    for j in (1,2,3):
        x = [Fraction(v[0]**2,v[j]**2) for v in coords]
        ranges.append((min(x),max(x)))
    floor = exclude_ports(ranges)
    ports = mixed_ports()
    time = q("26.81672869387187")
    errors = residuals(common_state(ports),time)
    assert all(v < q("3.18e-7") for v in errors)
    print("depth 16 convex-mixture exclusion above",floor.str(9))
    print("depth 17 fair mixture",counts,"33-copy common gauge residual bounds",
          *(v.upper().str(9) for v in errors))
    # With only the averaged ports constrained, mix the empty word with
    # the three positive norm-five generators. Each generator contributes
    # probability 4/5 to its designated nonscalar port.
    slots = 2**26
    axes = [round(Fraction(5*v*slots,4*sum(ports))) for v in ports[1:]]
    assert axes == [6681837,12419288,34240280] and sum(axes)<slots
    one_step = [5*slots-4*sum(axes), *(4*k for k in axes)]
    one_errors = residuals(common_state(one_step),time)
    assert all(v < q("3.18e-7") for v in one_errors)
    print("depth one with a bypass: maximum residual below",max(v.upper() for v in one_errors).str(9))


if __name__ == "__main__":
    verify()
