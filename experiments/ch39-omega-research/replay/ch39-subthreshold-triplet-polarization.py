"""Compute the finite momentum dependence of the coloured scalar below its threshold."""
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

matrix = run_path(str(Path(__file__).with_name("ch39-neutral-matrix-running.py")))
poles, geometry, q = (matrix[k] for k in ("poles","geometry","q"))
x = q(1).exp()
p = 2*matrix["flow"]["step"]


def scalar_integral(ratio):
    # r = q_E^2/M3^2. The coefficient is the exact beta integral of
    # (1-2u)^2 [u(1-u)]^n divided by n.
    assert 0 < ratio < 1
    total = arb(0)
    for n in range(1,25):
        term = ratio**n / (n*(2*n+3)*(2*n+1)*comb(2*n,n))
        total += term if n%2 else -term
    n = 25
    remainder = ratio**n / (n*(2*n+3)*(2*n+1)*comb(2*n,n))
    enclosure = total + arb(0).union(remainder)
    a = ratio/(4+ratio)
    closed = 2*(a.sqrt().atanh()/a.sqrt()-1)/(3*a)-q("2/9")
    return enclosure.intersection(closed)


def finite_trace(trace, energy, triplet):
    r = energy/triplet**2
    delta = 48*(scalar_integral(x*r)-scalar_integral(r))/(16*arb.pi()**2)
    assert delta > 0
    addition = arb_mat(3,3)
    addition[0,0] = delta/p
    return trace+addition


def dirac_difference(energy_ratio):
    """One colourless vectorlike Dirac fermion with unit family charge."""
    assert 0 < energy_ratio and x*energy_ratio < 1
    total = arb(0)
    for n in range(1, 15):
        moment = 1/q(n*(2*n+3)*comb(2*n+2,n+1))
        term = moment*(x**n-1)*energy_ratio**n
        total += term if n%2 else -term
    n = 15
    last = ((x**n-1)*energy_ratio**n)/(n*(2*n+3)*comb(2*n+2,n+1))
    return (total+arb(0).union(last))/(2*arb.pi()**2)


def verify():
    time = poles["study"]["nominal"][2]
    trace = matrix["charge_trace"]()
    targets = [poles["study"]["comp"]["model"]["target"](n,time) for n in (11430,11431,11432)]
    for case,logs in enumerate(poles["retune"]["log_sets"]):
        (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
        c1,c2,_,n1,_ = poles["reference_readings"](logs)
        e,b,d = n1**2/scale,c/rho,eta/rho
        n2 = geometry["second_neutral_mass"](scale,z,b,d)
        start = matrix["massive_exchange"](e,z,b,d)[2,2]
        P = matrix["massive_exchange"](e*x,z,b,d)
        gy = matrix["flow"]["running"]["hypercharge_window"](time+n1.log()).mid()
        T = P*(1/rho)+matrix["hypercharge"]*gy**2
        updated = finite_trace(trace,e,masses[0])
        ref = rho*(T.inv()-updated*p).inv()[2,2]/(x*start)
        old = rho*(T.inv()-trace*p).inv()[2,2]/(x*start)
        err = arb(0,q("1e-10").upper())
        data = [m*(1+err) for m in (c1,c2,n1,n2)]
        l,zz,_ = poles["charged_geometry"](*data[:2])
        bb,dd = geometry["neutral_shape"](l,zz,data[2:])
        ee = data[2]**2/l
        triplet,_ = geometry["physical_masses"](bb,arb(1),zz)
        denominator = x*matrix["massive_exchange"](ee,zz,bb,dd)[2,2]
        P = matrix["massive_exchange"](ee*x,zz,bb,dd)
        inverse_trace = finite_trace(trace,ee,triplet).inv()
        ks = [1/(poles["offset"]+t/(15*l)) for t in targets]
        low = matrix["enhancement"](P,inverse_trace,ks[0],q(1),p)/denominator
        high = matrix["enhancement"](P,inverse_trace,ks[2],q(0),p)/denominator
        old_max = matrix["enhancement"](P,trace.inv(),ks[1],q(1),p)/denominator
        print(case,"q/M3",(e*x).sqrt()/masses[0],"finite relative shift ppm",
              ((ref/old-1)*10**6).str(11),"distance from old upper edge ppm", ((ref/old_max-1)*10**6).lower().str(10))
        print("updated amplitude",ref.str(14),"outward lower margins ppm",
              ((1-low/ref)*10**6).lower().str(11), ((high/ref-1)*10**6).lower().str(11))
        assert q("27.44") < (ref/old-1)*10**6 < q("27.54")
        assert ref > old_max*(1+q("26.10e-6"))
        assert low < ref*(1-q("2.18e-6")) and high > ref*(1+q("1.43e-6"))
        # Use exact masses and the same hypercharge in the alternative model.
        # A vectorlike singlet's mass is free and it carries no SM charge.
        alternate_rho = poles["offset"]+targets[0]/(15*scale)
        alternate_T = matrix["massive_exchange"](e*x,z,b,d)*(1/alternate_rho)+matrix["hypercharge"]*gy**2
        projected = (alternate_T.inv()-updated*p).inv()[2,2]
        needed = 1/projected-alternate_rho/(ref*x*start)
        assert needed > 0
        # x e/M_3^2 < 1; the mass search starts at the triplet threshold.
        def gap(r):
            return dirac_difference(r)-needed
        r = poles["roots"](gap,q(".00001"),e/masses[0]**2)
        band = r+arb(0,q("1e-12").upper())
        assert gap(band.lower()) < 0 < gap(band.upper())
        mass_ratio = (e/band).sqrt()/masses[0]
        assert mass_ratio > 1 and abs(mass_ratio/mass_ratio.mid()-1) < q("1e-7")
        adjusted = alternate_rho*projected/(1-dirac_difference(band)*projected)/(x*start)
        assert adjusted.overlaps(ref)
        print("11430 vectorlike singlet mass / M3",mass_ratio.str(12),flush=True)



if __name__ == "__main__":
    verify()
