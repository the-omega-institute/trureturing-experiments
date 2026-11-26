"""Displacement and stiffness required when the singlet reference and matrix Higgs move together."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb

kin = run_path(str(Path(__file__).with_name("ch39-shell-gw-native-kinetic-cost.py")))
n = kin["n"]
m = kin["m"]
adj,count,rows = kin["adj"],kin["count"],kin["rows"]


if __name__=="__main__":
    t = n["relative_phase"](arb(676))
    y = kin["arc"]["conjugate_output"](t)
    D = adj(y)*y
    K = kin["scalar_kinetic"]()
    b,x,weight = arb(1020),arb(250),arb(16)
    edge = arb(".46")
    assert count(rows(K+b*D-x*m["G"]),edge)==0
    u2,d2 = m["u2"],m["d2"]
    # h bounds the kinetic quotient if the trial copy fixes the reference to q0.
    h = arb("1.347")
    slope = x-b*(1+weight)*d2/u2
    assert slope > 0
    excess = edge+slope*u2/(4*arb.pi())**2-h
    assert excess > 0
    gain = arb(".515")
    delta2 = excess/(b*(1+1/weight)*gain**2)
    assert delta2.sqrt() > arb(".030")
    # q0 is an exact ground state of the tangent Hamiltonian. Lift the native
    # reference spectrum to a bound on the energy cost of any change of its line.
    p,sing = m["c"]["p"],m["c"]["c"]["sing"]
    phi = -p[0].conjugate()*sing[0,3]*p[3]
    phases = [acb(0,(-1)**i*t/2).exp() for i in range(6)]
    corr = (phi*phases[1]/phases[0]).real
    beta = 676*(corr-arb(25)/26)
    pts = m["c"]["m"]["m"]["pts"]
    tangent = acb_mat([[p[i].conjugate()*sing[i,j]*p[j]*
                         (1-beta if sum(a*b for a,b in zip(pts[i],pts[j]))==25 else 1)
                         for j in range(6)] for i in range(6)])
    q0 = acb_mat([[z] for z in phases])
    e = m["hs"](q0,tangent*q0).real/6
    gap = arb("1.15")
    assert count(rows(tangent),e+gap)==1
    # The best global phase makes ||q-q0||² <= 2||(I-P0)q||².
    # Squaring the six unit phases expands each chord length by at most two.
    reference_cost = gap*delta2/8
    # Use local block gains and the best global phase of the squared reference.
    F = m["F"]
    block_gain = arb(".34")
    for i in range(6):
        site = acb_mat([[F[j,9*i+b] for b in range(9)] for j in range(9)])
        assert count(rows(site*adj(site)),block_gain**2)==9
    b,x,weight = arb(320),arb(110),arb(20)
    edge = arb(".8392")
    assert count(rows(K+b*D-x*m["G"]),edge)==0
    slope = x-b*(1+weight)*d2/u2
    assert slope>0 and edge+slope*u2/(4*arb.pi())**2>h
    stiffness = b*(1+1/weight)*4*block_gain**2/gap
    assert stiffness < arb("135.2")
    print("necessary squared-phase chord above",delta2.sqrt().lower().str(8),
          "excess reference energy above",reference_cost.lower().str(8),
          "necessary stiffness / matrix norm² below",stiffness.upper().str(8))
