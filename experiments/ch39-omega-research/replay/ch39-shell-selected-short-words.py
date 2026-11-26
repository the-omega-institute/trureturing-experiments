"""Use one or two uniform lifted shells to select norm-five preparation words."""
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement, product
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb

root = Path(__file__).parent
short = run_path(str(root/"ch39-short-word-gauge-mixture.py"))
local = run_path(str(root/"ch39-local-product-contact.py"))
gauge, shortest, model, q = (short[k] for k in ("gauge","shortest","model","q"))
events, words = short["events"],short["words"]
selected = ((527001,246580,336568,557790),(-527151,246378,336544,557752),
            (527281,246672,336138,557744),(-527051,246522,336154,558018),
            (-526967,246408,335994,558244),(527281,246672,336138,557744))


def radial_counts():
    sites, pole = local["lifted_sites"]()
    distance = local["geometry"]["squared"]
    levels = Counter(distance(v,pole) for v in sites)
    assert levels == {0:1,68:9,136:9,544:6}
    def group(site):
        dist = distance(site,pole)
        return 0 if dist==0 else 1 if dist==68 else 2
    pairs = Counter(tuple(sorted((group(a),group(b)))) for a,b in product(sites,repeat=2))
    return tuple(levels[d] for d in sorted(levels)), [pairs[key] for key in sorted(pairs)]


def scale_constraints(time):
    reference_mass = short["mixed_ports"]()
    ref = tuple(Fraction(v,sum(reference_mass)) for v in reference_mass)
    ref_state = short["common_state"](reference_mass)
    ys = [v.log() for v in ref_state[3:]]
    for _ in range(128):
        ref_state = shortest["reference"]["advance"](ref_state,time/128,14)
    target = model["target"](11431,time)
    errors = [1/(r*u)-target for r,u in zip(shortest["ratios"],ref_state)]
    mix = [a+b for a,b in events["mixing"]["displacement"](27)]
    volume = events["mixing"]["volume"]
    xl,xh,logs = volume["input_bounds"]()
    jac,_ = volume["variations"]([l.union(h) for l,h in zip(xl,xh)],logs,refinements=2)
    anchor = q("429/16")
    targets = model["target"](11431,anchor+q("1/8")).union(model["target"](11431,anchor-q("1/8")))
    speeds = short["floating"]["flow_bounds"]()
    beta = [1+(s-speeds[2])/(2*targets+speeds[2]) for s in speeds[:2]]
    tolerance = arb(2)**-21+q("1.3e-16")
    constraints = []
    for i,b in enumerate(beta):
        row = [jac[i,j]/shortest["ratios"][i]-b*jac[2,j]/shortest["ratios"][2] for j in range(12)]
        margin = mix[i]+b*mix[2]+(1+b)*tolerance
        constraints.append((row, errors[i]-b*errors[2], margin))
    return ref, ys, logs, constraints, beta, mix


def affine_port_bounds(ports, ref, ys, logs, constraints):
    # Linearize x_i = p_0/p_i. The remainder covers the convex hull.
    xs = [ref[0]/p for p in ref[1:]]
    unit = 5**17
    probs = [tuple(Fraction(a,unit) for a in row) for row in ports]
    terms = [[(v[0]-x*v[k])/ref[k] for k,x in enumerate(xs,1)] for v in probs]
    remainder = []
    for k in range(1,4):
        err = max(abs(v[k]-ref[k]) for v in probs)*max(abs(z[k-1]) for z in terms)/min(v[k] for v in probs)
        remainder.append(arb(0,q(err).upper()))
    bounds = []
    for row, center, margin in constraints:
        offset = center+arb(0,margin.upper())
        offset += sum(a*(l-v) for a,l,v in zip(row[3:],logs,ys,strict=True))
        offset += sum(a*rem for a,rem in zip(row[:3],remainder,strict=True))
        tests = []
        for affine in terms:
            test = 10**9*(offset+sum(a*q(z) for a,z in zip(row[:3],affine,strict=True)))
            tests.append((int(test.lower().floor().fmpz()),int(test.upper().ceil().fmpz())))
        bounds.append(tests)
    return np.array(bounds,dtype=np.int64)


def coarse_tests(ports, time):
    ref, ys, logs, constraints, beta, mix = scale_constraints(time)
    return affine_port_bounds(ports, ref, ys, logs, constraints), beta, mix


def verify():
    single, pair = radial_counts()
    library = [(v,d) for d in (16,17) for v,_ in words["hits"](d)]
    assert len(library)==43
    ports = [[x*x*5**(17-d) for x in v] for v,d in library]
    time = q("26.81672869387187")
    tests,beta,mix = coarse_tests(ports,time)
    assert int(np.max(np.abs(tests))) < 10**8
    same_size = list(combinations_with_replacement(range(43),2))
    middle = np.array(same_size)
    survivors = []
    assert single==(1,9,9,6)
    for i,j in product(range(43),repeat=2):
        out = tests[:,i:i+1,:]+6*tests[:,j:j+1,:]+9*(tests[:,middle[:,0],:]+tests[:,middle[:,1],:])
        mask = np.all((out[:,:,0]<=0)&(out[:,:,1]>=0),axis=0)
        survivors.extend((i,k,l,j) for k,l in middle[mask])
    rejected_by_table = 0
    excluded = 0
    for ids in survivors:
        mass = [sum(n*ports[j][axis] for n,j in zip(single,ids)) for axis in range(4)]
        values = short["table_readout"](mass)
        if any(not events["low"][k]<=v<=events["high"][k] for k,v in zip(events["names"],values)):
            rejected_by_table += 1
            continue
        state = [arb(v)/mass[0] for v in mass[1:]]+[q(v) for v in values]
        for _ in range(128):
            state = shortest["reference"]["advance"](state,time/128,14)
        target = model["target"](11431,time)
        residual = [1/(r*u)-target for r,u in zip(shortest["ratios"],state)]
        for i,b in enumerate(beta):
            if abs(residual[i]-b*residual[2]) > (1+b)*(arb(2)**-21+q("1.3e-16"))+mix[i]+b*mix[2]:
                excluded += 1
                break
        else:
            raise AssertionError(("unresolved single-shell assignment",ids))
    assert (len(survivors),rejected_by_table,excluded)==(3,2,1)
    print("one shell: all radial assignments excluded; refined",len(survivors),rejected_by_table,excluded)
    assert pair==[1,18,30,81,270,225] and all((v,17) in library for v in selected)
    mass = [sum(n*v[j]**2 for n,v in zip(pair,selected,strict=True)) for j in range(4)]
    trial = q("26.81668495751949")
    errors = short["residuals"](short["common_state"](mass),trial)
    assert all(v<q("3.18e-7") for v in errors)
    print("two shells: 625 equal pairs and all-mixing errors below",*(v.upper().str(9) for v in errors))


if __name__ == "__main__":
    verify()
