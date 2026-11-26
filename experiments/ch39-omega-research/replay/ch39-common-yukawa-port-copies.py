"""Certify a common-reference readout of nine Yukawa couplings on 33 copies."""
from collections import defaultdict
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

root = Path(__file__).parent
model = run_path(str(root/"ch39-shortest-matching-rotor.py"))
mixing = run_path(str(root/"ch39-ckm-gauge-displacement.py"))
exact = model["reference"]["rational"]
names = ("u","c","t","d","s","b","e","mu","tau")
centers = dict(zip(names,exact["YUKAWA"],strict=True))
low = {name:(c*Q(".99"))**2 for name,c in centers.items()}
high = {name:(c*Q("1.01"))**2 for name,c in centers.items()}


def orbit_weights(copies, coordinates=None):
    if coordinates is None:
        coordinates = model["match_coord"]
        assert sum(v*v for v in coordinates) == 5**21
    return orbit_weights_from_ports(copies, [v*v for v in coordinates])


def orbit_weights_from_ports(copies, ports):
    powers = [[v**k for k in range(copies+1)] for v in ports]
    fac = [factorial(k) for k in range(copies+1)]
    weights = defaultdict(int)
    for i in range(copies+1):
        for j in range(copies+1-i):
            for k in range(copies+1-i-j):
                scalar = copies-i-j-k
                mult = fac[copies]//(fac[scalar]*fac[i]*fac[j]*fac[k])
                key = tuple(sorted((i,j,k),reverse=True))
                weights[key] += mult*powers[0][scalar]*powers[1][i]*powers[2][j]*powers[3][k]
    assert sum(weights.values()) == sum(ports)**copies
    return weights


def possible_light_events(weights):
    # Work in integer mass units. All nine event masses plus the reference
    # must fit in the unit total; the three lightest bins must be disjoint.
    max_ref = sum(weights.values())/(1+sum(low.values()))
    light = ("e","u","d")
    pool = sorted(w for w in weights.values() if w<=max_ref*high["d"])

    def subsets(start,mass,mask,used,min_mass,max_mass):
        if mass and mass>=min_mass:
            yield mass,mask
        for j in range(start,len(pool)):
            if mass+pool[j]>max_mass:
                break
            if not used & (1<<j):
                yield from subsets(j+1,mass+pool[j],mask|(1<<j),used,min_mass,max_mass)

    def place(index,min_ref,max_ref,used):
        if index == len(light):
            return True
        name = light[index]
        for mass,mask in subsets(0,0,0,used,min_ref*low[name],max_ref*high[name]):
            a = max(min_ref,mass/high[name])
            b = min(max_ref,mass/low[name])
            if a<=b and place(index+1,a,b,used|mask):
                return True
        return False

    return place(0,0,max_ref,0)


def event_masses(weights):
    masses = defaultdict(int)
    used = set()
    table = (root/"ch39-common-yukawa-events.tsv").read_text().splitlines()
    assert table[0].split() == ["largest","middle","smallest","event"]
    for row in table[1:]:
        a,b,c,name = row.split()
        key = (int(a),int(b),int(c))
        assert key not in used
        used.add(key)
        masses[name] += weights[key]
    assert set(masses) == set(names)|{"ref"}
    return masses,used


def verify():
    ctx.prec = 128
    for copies in range(1,33):
        assert not possible_light_events(orbit_weights(copies)), copies
    weights = orbit_weights(33)
    masses,used = event_masses(weights)
    squares = [Q(masses[n],masses["ref"]) for n in names]
    assert all(low[n]<=z<=high[n] for n,z in zip(names,squares))
    a,*imag = model["match_coord"]
    state = [(arb(v)/a)**2 for v in imag]+[model["q"](z) for z in squares]
    time = model["q"]("26.8162148027")
    for _ in range(128):
        state = model["reference"]["advance"](state,time/128,14)
    change = [a+b for a,b in mixing["displacement"](27)]
    target = model["model"]["target"](model["match_label"],time)
    errors = [abs(1/(u*r)-target)+d for u,r,d in zip(state,model["ratios"],change)]
    assert all(e < model["tolerance"] for e in errors)
    # Optional discarded orbits can be assigned to the bottom channel. Bound
    # their total mass and propagate the resulting log-Yukawa change through
    # the diagonal flow before using the same uniform mixing bound.
    spare = sorted(w for key,w in weights.items() if key not in used)
    extra = sum(spare[:139])
    assert extra*200 < masses["b"]
    assert Q(masses["b"]+extra,masses["ref"]) <= high["b"]
    v = model["volume"]
    xl,xh,logs = v["input_bounds"]()
    jac,_ = v["variations"]([l.union(h) for l,h in zip(xl,xh)],logs,refinements=2)
    shift = (1+model["q"](Q(extra,masses["b"]))).log()
    family = [err+abs(jac[i,8])*shift/r for i,(err,r) in enumerate(zip(errors,model["ratios"]))]
    assert all(err < model["tolerance"] for err in family)
    print("Optional bottom-channel assignments: 2^139, family residual upper bound",
          max(err.upper() for err in family).str(10))
    print("33-copy Yukawa factors",*[round(float(z/c**2)**.5,10) for z,c in zip(squares,centers.values())])
    print("Minimal common-reference copies 33, gauge residual for all mixing <",
          max(e.upper() for e in errors).str(10))


if __name__ == "__main__":
    verify()
