"""Continue the neutral logarithms through the triplet threshold up to the first W_R."""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
running = run_path(str(root/"ch39-neutral-abelian-running.py"))
flavour, study, q = running["flavour"], running["study"], running["q"]
retuned, comp = study["retuned"], study["comp"]
adj = study["adjoint"]


def graph_gap():
    neutral = retuned["current"]["neutral"]
    def checked(*args, **kw):
        return neutral["pole_polynomials"](*args,**kw,gap_denominator=8)
    unit,_,_,_ = neutral["graph_constants"](port_builder=checked)
    layers = ((1,5,29),(7,17,23),(11,11,25),(17,17,17))
    vertices = [p for p in product(range(30),repeat=3)
                if sum(v*v for v in p)==867 and tuple(sorted(p)) in layers]
    def left_edge(p,r):
        return neutral["graph"]["squared"](p,r)!=136 or (11,11,25) in (tuple(sorted(p)),tuple(sorted(r)))
    checked(vertices,(32,152,56,136),27,edge_filter=left_edge)
    return unit/8


def right_threshold(masses,c,z,rho,eta,gap):
    triplet,octet = masses
    unit,colour,right = study["ports"][-3:]
    # g^2 w^2/2=1/(4z) and 16g^2 v_Sigma^2=d^2/h^2.
    shift_r = 1/(4*z)
    h2 = triplet**2*18*rho/c
    shift_c = shift_r+(octet**2-triplet**2)/h2
    root_r = study["enclose_roots"](lambda e:1+shift_r*right(e,unit),".02",".035")
    root_c = study["enclose_roots"](lambda e:1+shift_c*colour(e,unit),".035",".05")
    assert triplet**2<root_r<root_c<gap
    assert min(c/eta,octet**2,octet**2-triplet**2)>root_r
    return root_r.sqrt()


def window_data(model,gap):
    L,z,rho,eta,c,masses = model
    cal = study["nominal"][2]+L.log()/2
    neutral = run_path(str(root/"ch39-neutral-current-threshold-inverse.py"))
    pole,charges,_,_ = neutral["first_window"](*masses,c,z/rho,eta,1/rho)
    t0 = cal+pole.log()/2
    h3,h8 = [cal+m.log() for m in masses]
    right = right_threshold(masses,c,z,rho,eta,gap)
    end = cal+right.log()
    gY = running["hypercharge_window"](t0)
    first = running["run_window"](h3-t0,gY,charges)
    yy = q("41/6")
    gY3 = gY/(1-(h3-t0)*yy*gY**2/(8*arb.pi()**2)).sqrt()
    second = running["run_window"](end-h3,gY3,first[4],triplet=True)
    print("W_R / M3", (right/masses[0]).str(12))
    return t0,h3,end,h8,first,second


def bounds(index,upper,windows,time,*,integrals=None):
    onset,h3,end,h8,first,second = windows
    state = flavour["before"](index,upper)
    pick = 2 if upper else 3
    # Split at each uncertain boundary; the gauge and flavour effects need
    # their own monotonicity bounds in intervals where the field can switch.
    active = [[arb(0)]*3,first[pick],second[pick],[arb(0)]*3]
    edges = [q(23),onset.lower(),onset.upper(),h3.lower(),h3.upper(),end.lower(),end.upper(),
             h8.lower(),h8.upper(),time]
    for a,b in zip(edges,edges[1:]):
        if a == b:
            continue
        assert a < b
        mid = (a+b)/2
        possible = [v for i,v in enumerate(active) if
                    (i==0 or mid>(onset,h3,end)[i-1].lower()) and
                    (i==3 or mid<(onset,h3,end)[i].upper())]
        rate = [min(v[i].lower() for v in possible) if upper
                else max(v[i].upper() for v in possible) for i in range(3)]
        scalar = [mid>(x.upper() if upper else x.lower()) for x in (h3,h8)]
        flow = {**adj["flow"](*scalar),"NEUTRAL_YUKAWA":rate}
        n = 1+int(12*float((b-a).upper()))
        for _ in range(n):
            state,tube,*values = comp["advance"](state,(b-a)/n,upper,flow=flow,need_tube=True,
                                                  integrate_state=integrals is not None)
            if integrals is not None:
                integrals.append(((b-a)/n,flow,values[0]))
            assert all(abs(v)<2 for v in tube)
    return [1/(r*v)-f-s for r,v,f,s in zip(adj["ratios"],state,first[0],second[0])]


def verify_speeds(windows):
    gauge = arb(0).union(arb(2))
    max_rate = max(*windows[-1][1],*windows[-2][1])
    for flow in (adj["flow"](False,False),adj["flow"](True,False),adj["flow"](True,True)):
        speed = [-2*(beta+sum(b*gauge-d*3*gauge for b,d in zip(bs,ds))/comp["loop"])
                 /(r*comp["loop"]) for beta,bs,ds,r in
                 zip(flow["BETA"],flow["B"],flow["D"],adj["ratios"],strict=True)]
        assert speed[1]>max_rate and speed[1]-speed[0]>max_rate and speed[2]-speed[1]>max_rate


def verify():
    fit = run_path(str(root/"ch39-finite-current-width-inverse.py"))
    _,_,L,z,rho,eta,c,masses = fit["inferred_parameters"]()
    windows = window_data((L,z,rho,eta,c,masses),graph_gap())
    print("second-window gauge shifts",*(v.str(10) for v in windows[-1][0]))
    cal = study["nominal"][2]
    verify_speeds(windows)
    for index in range(3):
        at = cal-q(".00008")
        lo,hi = [bounds(index,upper,windows,at) for upper in (False,True)]
        radius = min(lo[1]-hi[0],lo[1]-hi[2])/2
        assert radius > q(("1.1813e-5","1.1343e-5","8.3708e-6")[index])
        for upper in (False,True):
            bounds(index,upper,windows,q(33))
        print("all-scale distance to any common value, group",index,radius.lower().str(10))


if __name__ == "__main__":
    verify()
