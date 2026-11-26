"""Global singlet phase vacuum with native Berry transport and a shortest-edge comparison."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb

arc = run_path(str(Path(__file__).with_name("ch39-shell-gw-arc-reference.py")))
m = arc["m"]
c = m["c"]["c"]
p,sing = m["c"]["p"],c["sing"]
hs,adj,count,rows = m["hs"],m["adj"],m["count"],m["rows"]


def relative_phase(kappa):
    def energy(signs):
        q = acb_mat([[p[i]*signs[i]] for i in range(6)])
        return hs(q,sing*q).real/6
    even = energy([1]*6)
    odd = energy([(-1)**i for i in range(6)])
    B,A = (odd-even)/2,(even+odd)/2
    phase = -p[0].conjugate()*sing[0,3]*p[3]
    c0 = arb(25)/26
    assert arb("-.0780") < phase.arg() < arb("-.0779")
    bmax = B*phase.real-B*abs(phase.imag)*c0/(1-c0**2).sqrt()
    assert 0 < bmax < arb("1.445")
    # In the two complementary C3 sectors the first native eigenvalue exceeds 3.71.
    # The other invariant-orbit eigenvalue is odd > 4 and the ground state is even.
    assert even < arb("3.71") < odd and count(rows(sing),arb("3.71")) == 1
    assert arb("3.71")-bmax > A
    def gap(b):
        d = B-b*phase
        return kappa*((phase*d.conjugate()).real/abs(d)-c0)-b
    lo,hi = arb(0),arb("1.445")
    assert gap(lo)>0 and gap(hi)<0
    for _ in range(65):
        mid = (lo+hi)/2
        if gap(mid)>0:
            lo = mid
        else:
            assert gap(mid)<0
            hi = mid
    return (B-lo.union(hi)*phase).arg()


if __name__=="__main__":
    t = relative_phase(arb(26)**2)
    assert arb(".1925047118") < t < arb(".1925047120")
    target = arc["central_target"]()
    paired = acb_mat(rows(m["F"])+rows(arc["conjugate_output"](t)))
    field = adj(paired)*(paired*adj(paired)).solve(target)
    strength = hs(field,field).real.sqrt()
    assert strength < arb("11.6") < 4*arb.pi()
    print("native kinetic plus normalized shortest-edge comparison: unique phase",
          t.str(11),"equal-coupling cost",strength.str(10))
