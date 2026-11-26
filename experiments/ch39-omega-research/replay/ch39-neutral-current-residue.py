"""Isolate the family-current residue and its two-loop leading logarithm."""
from pathlib import Path
from runpy import run_path

from flint import arb

neutral = run_path(str(Path(__file__).with_name("ch39-shared-quartic-neutral.py")))


def verify():
    g,u,w,m3,_,active = neutral["parameters"]()
    c = 18*u*u
    def positive(energy):
        return -1/(2*w*w)-active(energy)
    limit = neutral["bisect"](lambda e:4/c-positive(e),".0001",".003")
    fraction_floor = 4/(c*limit*active(limit,True))
    def measure(energy):
        h = positive(energy)
        return c*h*h/(c*h*h+4*(energy*active(energy,True)-h))
    assert (measure(limit)-fraction_floor).contains(0)
    def full(energy):
        return active(energy)+1/(2*w*w)+4*g*g/(c*g*g-energy)
    assert c*g*g>limit
    energy = neutral["bisect"](full,".0001",limit)
    fraction=measure(energy)
    assert fraction_floor<fraction<1
    g_eff_sq=energy*fraction/c
    assert (energy*positive(energy)/(c*positive(energy)-4)-g*g).contains(0)
    log_range=(m3/energy.sqrt()).log()
    shift=48*g_eff_sq*log_range/(16*arb.pi()**2)**2
    assert arb(".95043290") < fraction_floor < arb(".95043291")
    assert arb(".00004445174") < shift < arb(".00004445175")
    tolerance=arb(2)**-21
    prefactor=48/(c*(16*arb.pi()**2)**2)
    def logarithm(e):
        return (m3/e.sqrt()).log()
    assert limit*arb(1).exp()<m3*m3
    # E^2 A' increases and E H decreases by spectral positivity. The
    # graph gap is above m3^2, as certified in the mass calculation.
    for scale_error,cut,width,copies in [(arb(1),arb(".000004"),".008067559276",5481),
                                        (arb(5)/96,arb(".00015"),".049886978042",144)]:
        assert scale_error*prefactor*fraction_floor*cut*logarithm(cut)>tolerance
        loss=(1-fraction_floor)*(2*cut/(m3*m3-cut)+2*cut*active(cut,True)/positive(cut))
        assert loss < arb(".12") < 1-1/(2*logarithm(cut))
        boundary=neutral["bisect"](lambda e:scale_error*prefactor*measure(e)*e*logarithm(e)-tolerance,
                                    ".000001",cut)
        h=positive(boundary)
        coupling_cut=(boundary*h/(c*h-4)).sqrt()
        assert abs(coupling_cut-arb(width))<arb("1e-12")
        assert copies-1 < g*g/(12*coupling_cut**2) < copies
        print("family-coupling cutoff",coupling_cut.str(12),"diagonal SU(3) copies",copies)
    print("all couplings: light-pole contact fraction infimum",fraction_floor.str(13))
    print("equal site and family couplings: mass / Q",energy.sqrt().str(12),
          "contact fraction",fraction.str(12),"family coupling",g_eff_sq.sqrt().str(12),
          "two-loop family-only leading log to first triplet",shift.str(12))


if __name__=="__main__":
    verify()
