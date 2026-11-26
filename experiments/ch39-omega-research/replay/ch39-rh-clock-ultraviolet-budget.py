"""Minimum clock unit for a matching with all alpha_i below one through a given energy."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

counts = run_path(str(Path(__file__).with_name("ch39-rh-clock-solution-count.py")))
fixed = counts["fixed"]
positive = fixed["positive"]
H = 125**2


def static_residual(w, data):
    _, _, p, modes, k, f = data
    q, mf, u, v = fixed["scales"](w, k)
    e = sum(n*(q/(H+lam*p*mf**2).sqrt()).log() for lam, n in modes)
    return 22*e-402*u-46*v-15*f


def matched_w(epsilon, data):
    # For the light-only branch G_static <= G <= G_static + 22 epsilon.
    def isolate(shift):
        lo, hi = arb("8.19"), arb("8.20")
        assert static_residual(lo, data)+shift < 0 and static_residual(hi, data)+shift > 0
        for _ in range(65):
            mid = (lo+hi)/2
            value = static_residual(mid, data)+shift
            if value > 0:
                hi = mid
            elif value < 0:
                lo = mid
            else:
                break
        return lo.union(hi)
    return isolate(0).union(isolate(22*epsilon))


def threshold_interval(exposure, budget, smallest, energy):
    lo, hi = int(smallest.unique_fmpz()), int(energy.unique_fmpz())
    for _ in range(68):
        if hi-lo <= 1:
            break
        mid = (lo+hi)//2
        val = exposure(mid)-budget
        if val > 0:
            lo = mid
        elif val < 0:
            hi = mid
        else:
            break
    assert exposure(lo) > budget and exposure(hi) < budget
    return lo, hi


def main():
    ctx.prec, ctx.cap = 1536, 3
    data = fixed["inputs"]()
    _, _, p, modes, k, f = data
    d1, d2, d3 = [8*arb.pi()**2/arb(g)**2 for g in (".467453", ".63811", "1.0583")]
    for energy in (arb("1e17"), arb("1.2209e19")):
        smallest = energy/1000
        assert smallest > arb("2.75e12")
        poles, rho, *_ = counts["at_unit"](data, smallest)
        t = (energy/1000).log()
        allowance = 6*(d2-2*arb.pi())+19*t
        def pure_load(unit):
            return 48*sum(positive((energy*smallest/(unit*c.sqrt())).log()) for c in poles)
        assert pure_load(smallest) > allowance
        qhi, mfhi, _, _ = fixed["scales"](arb("8.19"), k)
        qlo, mflo, _, _ = fixed["scales"](arb("8.20"), k)
        c1 = poles[0]
        lmax = 6*mfhi**2
        xmax = H+p*lmax
        assert xmax < qlo**2 and qhi**2 < c1
        ratio = lmax/(c1-xmax)
        epsilon = arb(23)/2*(1+ratio).log()
        assert (H+min(lam for lam, n in modes)*p*mflo**2)/(1+ratio) > 1000**2
        w = matched_w(epsilon, data)
        q, mf, u, v = fixed["scales"](w, k)
        light_load = (15*f+402*u+46*v)/22+23*(energy/q).log()
        # Heavy shifted roots interlace above the poles. Their total logarithmic
        # displacement is fixed by the determinant and bounded by epsilon.
        def load(unit):
            return light_load+pure_load(arb(unit))-arb(0).union(epsilon)
        # The 23 light roots have very small reference component at the crossing.
        ref = lmax/c1/(1-xmax/c1)**2
        emin = (1-arb(1)/64)/(1+ref)
        rw = 22*23*((1+arb(4)/11)*emin-arb(4)/11)-448
        assert rw > 0
        drift = 22*23*ref/rw
        assert arb(15)/11*23*drift < 25
        # Above the last possible pure-state exposure the whole weak load is safe.
        assert light_load < allowance and q < energy
        lo, hi = threshold_interval(load, allowance, smallest, energy)
        assert arb(hi)*c1.sqrt()/smallest < energy
        # Once the weak-channel bound is met, it also bounds hypercharge.
        trip = positive((energy/(mflo*arb(8).sqrt())).log())
        sing = positive((energy/(mflo*arb(24).sqrt())).log())
        hyper_floor = d1-arb(41)/10*t-allowance/10-arb(16)/5*trip-arb(4)/5*sing
        print("benchmark", energy.str(8), "weak-unit bracket", lo, hi,
              "hypercharge floor on weak boundary", hyper_floor.str(9), flush=True)
        # Hypercharge can be the tighter condition. Its load budget is known
        # from the same w interval, and its inverse increases while a pure pole
        # is open even after allowing for motion of the fermion thresholds.
        trip_box = positive((energy/(mf*arb(8).sqrt())).log())
        sing_box = positive((energy/(mf*arb(24).sqrt())).log())
        hyper_allow = 10*(d1-2*arb.pi())-41*t-32*trip_box-8*sing_box
        assert arb(15)/11*(23+40)*drift < 25 and light_load < hyper_allow
        lo1, hi1 = threshold_interval(load, hyper_allow, smallest, energy)
        assert arb(hi1)*c1.sqrt()/smallest < energy
        print("hypercharge-unit bracket", lo1, hi1, "overall", max(lo,lo1), max(hi,hi1))
        assert hi < lo1
        assert (d3+7*(q/1000).log()-2*u) > 2*arb.pi()


if __name__ == "__main__":
    main()
