"""A light-heavy mixing that keeps the 30 heavy masses gapped on a finite disk."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpq, fmpq_mat

ROOT = Path(__file__).parent


def check_minors():
    # Interpolation on this grid proves the identity in z,s for real z,s
    # (degree <=6 in each). For complex z,s the same minors are monomials.
    for z in range(-3, 4):
        for s in range(-3, 4):
            mass = fmpq_mat([[z, 0, s, 0], [s, z, 0, 1], [0, s, z, 0]])
            minors = [fmpq_mat([[mass[i, j] for j in range(4) if j != drop]
                                  for i in range(3)]).det() for drop in range(4)]
            assert minors == [s*s, -z*z, -s*z, z**3+s**3]
            assert (mass*mass.transpose()).det() == sum(c*c for c in minors)
    # r+s = 25/24 - (2/3)(r-3/4)^2 for s=(2/3)(1-r^2).
    for r in [fmpq(i, 4) for i in range(5)]:
        assert r+fmpq(2, 3)*(1-r*r) == fmpq(25, 24)-fmpq(2, 3)*(r-fmpq(3, 4))**2


def check_selected_gap():
    ext = run_path(str(ROOT/"ch39-selected-extended-flow.py"))
    v, state = ext["initial"]()
    old_x, old_y = [v*t.sqrt() for t in state[2:4]]
    mu = old_y/arb(2).root(6)
    assert (2*mu**6-old_y**6).contains(0)
    # The root is m_x'/m_x; the comparable endpoints certify the small shift.
    shift = (mu/old_x)**2
    def residual(t):
        w = t*t
        return w*w*(w+shift)-1
    lo = 1-arb("2.82e-10")
    hi = 1-arb("2.80e-10")
    assert residual(lo)<0<residual(hi)
    gap = mu/(4*(1+(arb(25)/24)**2))
    assert arb("2.7263e10")<mu<arb("2.7264e10")
    assert arb("3.2688e9")<gap<arb("3.2690e9")
    print("Mass mixing / GeV",mu.str(12),"uniform heavy singular-value bound / GeV",gap.str(12))


def check_transport():
    for count in (3, 6, 9, 19):
        angle = arb.pi()/(10*count)
        link = (1+acb(0, -2*angle).exp())/2
        product = link**count
        answer = acb(0, -arb.pi()/10).exp()
        assert (product-answer*angle.cos()**count).contains(0)
        assert 0 < abs(product) and abs(product) < 1
    print("Equal mixing on the Q coordinates gives projected loop phase exp(-i*pi/10).")


if __name__ == "__main__":
    ctx.prec = 160
    check_minors()
    check_selected_gap()
    check_transport()
