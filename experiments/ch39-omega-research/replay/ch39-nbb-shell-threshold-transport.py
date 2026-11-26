"""Transport the joint, signed-band filter to the whole nine-point shell and the one-loop crossing."""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

directory = Path(__file__).parent
joint = run_path(str(directory/"ch39-rh-shell-clock-step.py"))
coherent = run_path(str(directory/"ch39-finite-nbb-coherent-four-windows.py"))


def mass_ratio(levels):
    def h(i,j):
        return (i*levels[i]+j*levels[j])/((i+j)*(levels[i]*levels[j]).sqrt())
    a,b,c = h(0,1),h(1,2),h(0,2)
    det = 1+2*a*b*c-a*a-b*b-c*c
    assert det>0 and c<a and c<b
    return ((1-c*c)/det).sqrt()


def main():
    ctx.prec, ctx.cap = 480,3
    levels,_,heights,_,weights = joint["joint_data"]()
    heights = [heights[j-1] for j in joint["support"]]
    error,_ = coherent["bound_filter"](heights,weights)
    # The four shell cosines have Lipschitz constants <= 2 log(6)<4.
    assert 2*arb(6).log()<4
    offdiag = 4*error
    # For each shell vector x_u=(1,u)/sqrt(2), ||x_u||_1<=2.
    # Its Gram perturbation in a normalized pair is <= 4*offdiag
    # before normalization, and <= 2*that/(1-that) afterwards.
    entry_error = 4*offdiag
    shell_error = 2*entry_error/(1-entry_error)
    assert shell_error < arb(2)**-149
    points = [v for v in product(range(5),repeat=3) if sum(x*x for x in v)==17]
    table = {(17+sum(x*y for x,y in zip(u,v))) for u in points for v in points if u!=v}
    assert len(points)==9 and 28 in table and max(table)==33
    # Use the rational record bins also for the reference-field spectrum.
    widened = [levels[0]]+[v+arb(0,arb(2)**-155) for v in levels[1:4]]
    for j in range(1,4):
        for i in range(j):
            h=(i*widened[i]+j*widened[j])/((i+j)*(widened[i]*widened[j]).sqrt())
            rec=coherent["record"]
            assert (h*rec["den"]).floor().unique_fmpz()==rec["bins"][j][i]
    base,filtered=mass_ratio(levels),mass_ratio(widened)
    assert abs(filtered-base)<arb("5.2e-36") and abs((filtered/base).log())<arb("1e-38")
    assert filtered>arb("528.89869680") and filtered<arb("528.89869682")
    # Differentiate the equations in ch39-finite-rh-pivot-calibration.py
    # in log(M_S) on its certified central bracket. Z'=-sum(h/(h+lambda*M_S^2)).
    eig = [arb(6)]*7+[arb(4)]*6
    for k,n in ((5,3),(3,2)):
        eig += [3-arb(k).sqrt()]*n+[3+arb(k).sqrt()]*n
    scalar = arb("2034928.14").union(arb("2034928.16"))
    positive = sum(125**2/(125**2+v*scalar**2) for v in eig)
    r_slope=(2552-1320*positive)/27808
    q_slope=1-(19712+352*positive)/27808
    assert r_slope>arb(".09") and q_slope>0 and q_slope<arb(".30")
    assert arb("3.9e12")*(arb("1e-38")*q_slope/r_slope).expm1()<arb("1.5e-25")
    print("Nine-shell Gram entries differ by <2^-149 from their rational targets; shortest arc 33/34, marked arc 14/17.")
    print("Change in M_F/M_S <5.2e-36, static one-loop Q shift <1.5e-25 GeV at the same central inputs.")
    print("The propagating-reference crossing stays fixed within the certified closed-fermion branch.")


if __name__ == "__main__":
    main()
