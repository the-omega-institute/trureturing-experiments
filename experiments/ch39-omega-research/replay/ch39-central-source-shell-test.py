"""Compare the shell ground sectors with the book's central low-energy input."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

leads = run_path(str(Path(__file__).with_name("ch39-unit-lead-gauge-shell.py")))
floating = leads["floating"]
gauge = floating["gauge"]


def central_landscape():
    speeds = floating["flow_bounds"]()
    shortest = gauge["shortest"]
    model = shortest["model"]
    t0 = arb("26.81676714125598")
    values = [x*x for x in model["centers"]+model["yukawa_centers"]]
    for _ in range(128):
        values = shortest["reference"]["advance"](values,t0/128,14)
    shifts = [a+b for a,b in gauge["events"]["mixing"]["displacement"](27)]
    initial = [1/(r*u)+arb(0,e.upper()) for r,u,e in zip(shortest["ratios"],values[:3],shifts,strict=True)]
    assert speeds[0] < speeds[1] < speeds[2]
    crossing = (initial[0]-initial[2])/(speeds[2]-speeds[0])
    mismatch = (initial[0]-initial[1]+crossing*(speeds[0]-speeds[1]))/2
    assert 23 < t0+crossing < 33 and arb(".00285") < mismatch < arb(".00293")
    # On each side of the crossing of X1 and X3, respectively X1-X2 and
    # X3-X2 bound the spread and grow away from the crossing.

    def penalty(shift, target, common, contrast):
        xs = [v+d*shift for v,d in zip(initial,speeds)]
        mean = sum(xs)/3
        spread = sum((xs[i]-xs[j])**2 for i,j in ((0,1),(0,2),(1,2)))/3
        return 3*common*(mean-target*(-2*shift).exp())**2+contrast*spread

    def minimum(label, common, contrast):
        target = model["target"](label,t0)
        mean = sum(initial)/3
        trial_shift = ((target/mean).log()/2).mid()
        shifted_mean = mean+sum(speeds)*trial_shift/3
        trial_shift = ((target/shifted_mean).log()/2).mid()
        upper = penalty(trial_shift, target, common, contrast).upper()
        # Only slabs whose lower bounds could compete with the upper bound
        # are refined. The exact minimizer is somewhere in this finite cover.
        bounds = []
        pending = [(arb(23)-t0,arb(33)-t0)]
        while pending:
            lo,hi = pending.pop()
            box = lo.union(hi)
            cost = penalty(box,target,common,contrast)
            if cost > upper:
                continue
            if hi-lo < arb("1e-5"):
                bounds.append(cost.lower())
                continue
            middle = (lo+hi).mid()/2
            pending.extend(((lo,middle),(middle,hi)))
        return min(bounds).union(upper)

    return minimum, initial, t0, mismatch


def verify():
    minimum, _, _, mismatch = central_landscape()
    pen = floating["penalty"]
    score = pen["competition"]["score"]
    free = pen["doublet"]["joint"]["contact"]["clock"]["free_ground"]
    bare = free(score(10762),2)-free(score(11431),2)
    assert arb("8.11e-6") < bare < arb("8.13e-6")
    # The opposite direction cannot beat the forward score in either label.
    p = pen["competition"]["pulses"]
    for n in (10762,11431):
        z = acb.zeta_zero(n)
        a,b = p["port_outputs"](p["generator"](z.imag),2,impulse=True)
        assert abs((-7*a[0]+24*b[0])/25)**2 > abs((-7*a[1]+24*b[1])/25)**2
    for count in (1,2,3):
        common,contrast = map(arb,leads["terminal_modes"](count))
        delta = bare+minimum(10762,common,contrast)-minimum(11431,common,contrast)
        if count==1:
            assert delta > arb("2e-6")
        elif count==2:
            assert delta < arb("-1.5e-6")
        else:
            assert delta < arb("-5.5e-6")
        print("leads",count,"E10762-E11431",delta.lower().str(8),delta.upper().str(8))
    print("central inputs: best possible common-target max error",
          mismatch.lower().str(9),mismatch.upper().str(9))


if __name__ == "__main__":
    verify()
