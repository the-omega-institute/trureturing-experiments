"""A cyclic singlet phase reference from the shortest arc and its equal-coupling cost."""
from pathlib import Path
from runpy import run_path
from itertools import combinations

from flint import acb, acb_mat, arb

m = run_path(str(Path(__file__).with_name("ch39-shell-gw-conjugate-cost.py")))
F,G = m["F"],m["G"]
adj,hs,rows,count = (m[x] for x in ("adj","hs","rows","count"))
a = arb
signs = [(-1)**(i//9) for i in range(54)]


def conjugate_output(t):
    phases = {e:acb(0,-e*t).exp() for e in (1,-1)}
    return acb_mat([[F[i,j].conjugate()*phases[e] for j,e in enumerate(signs)] for i in range(9)])


def central_target():
    ref = m["m"]["orb"]["reference"]()
    ck,du,dd = ref
    import numpy as np
    U,A,B = [m["m"]["four"]["exact_frame"](np.array(z).reshape(3,3)) for z in m["units"]]
    u,d = U*du*A*m["up"][-1],U*ck*dd*B*m["down"][-1]
    return acb_mat([[z] for z in u.entries()+d.conjugate().entries()])


def tiny_twist_exclusion():
    # D(G-62 D* conj(G) D)D* = D G D* - 62 conj(G), with ||D'|| <= 1.
    dimension = len(signs)
    bound = (m["u2"]-62*m["d2"])/(4*a.pi())**2
    step = a(".02")
    assert count(rows(G),a(".27")) == dimension
    # Linear interpolation of a Hermitian curve with ||P''|| <= 4||G||.
    limit = bound-4*a(".27")*step**2/8
    for j in range(-6,7):
        t = j*step
        e = acb(0,2*t).exp()
        rot = acb_mat([[G[i,k]*(1 if si==sk else e if si<sk else e.conjugate())
                        for k,sk in enumerate(signs)] for i,si in enumerate(signs)])
        assert count(rows(rot-62*G.conjugate()),limit)==dimension


def arc_reference():
    pts = m["c"]["m"]["m"]["pts"]
    inner = {(i,j):sum(u*v for u,v in zip(pts[i],pts[j])) for i,j in combinations(range(6),2)}
    assert max(inner.values()) == 25
    assert {edge for edge,z in inner.items() if z==19} == {(i,j) for i,j in inner if i%2==j%2}
    assert sorted(v for edge,z in inner.items() if z==25 for v in edge)==list(range(6))
    t = (a(25)/26).acos()
    target = central_target()
    for sign in (1,-1):
        paired = acb_mat(rows(F)+rows(conjugate_output(sign*t)))
        v = adj(paired)*(paired*adj(paired)).solve(target)
        cost = hs(v,v).real.sqrt()
        assert cost < a("9.82") < 4*a.pi()
        print("shortest arc sign",sign,"equal-coupling strength",cost.str(10))
    sing = m["c"]["c"]["sing"]
    p = m["c"]["p"]
    def energy(vector):
        v = acb_mat([[z] for z in vector])
        return hs(v,sing*v).real/6
    delta = (energy([(-1)**i*z for i,z in enumerate(p)])-energy(p))*(t/2).sin()**2
    assert a(".07668") < delta < a(".07670")
    print("reference kinetic cost per site above ground",delta.str(9))


if __name__ == "__main__":
    tiny_twist_exclusion()
    arc_reference()
