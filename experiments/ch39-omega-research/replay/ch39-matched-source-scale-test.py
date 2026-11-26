"""Reduce the full matched-ground test to stationarity in the 11431 sector."""
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
rigidity = run_path(str(root/"ch39-matched-ground-rigidity.py"))
leads, floating = rigidity["leads"], rigidity["floating"]


def verify():
    # The expanded 10400--12499 clock comparison is independent of its gauge source.
    leads["verify"]()
    target = floating["gauge"]["shortest"]["model"]["target"]
    lead_modes = [tuple(map(arb,leads["terminal_modes"](k))) for k in (1,2,3)]
    speeds = floating["flow_bounds"]()
    tolerance = arb(2)**-21
    error = arb("1.3e-16")
    anchor = arb(429)/16
    band = anchor+arb(0,arb(1)/8)
    # Reuse one exponential after taking differences of the boundary targets.
    unit = arb(27)
    factor = (-2*(band-unit)).exp()
    radius = arb(".061")
    low_difference = (target(10399,unit+radius)-target(11431,unit))*factor
    high_difference = (target(11431,unit)-target(12500,unit-radius))*factor
    distance = min(low_difference.lower(),high_difference.lower())-2*error-tolerance-arb(".09")*radius
    spread = arb(".1435")*radius-2*tolerance
    mean_rate = sum(speeds)/3
    centered = [abs(v-mean_rate).upper() for v in speeds]
    low = target(11431,anchor+arb(1)/8)-error
    high = target(11431,anchor-arb(1)/8)+error
    assert distance > 0 and spread > 0 and low > arb("2.58")
    for count,(m,d) in enumerate(lead_modes,1):
        trial = 3*d*tolerance**2
        exterior = min((3*m*distance**2).lower(),(d*spread**2/2).lower())-arb(".000011193")-trial
        assert exterior > arb("3.29e-6") and arb("8e-6")-trial > arb("7.99e-6")
        ratio = d/m
        narrow = (2*(3*ratio).sqrt()+2)*tolerance/arb(".1435")
        qlo,qhi = low*(-2*narrow).exp(),high*(2*narrow).exp()
        slope = 2*qlo-arb(".09")
        residual = tolerance+arb(".09")*narrow+high*((2*narrow).exp()-(-2*narrow).exp())
        curvature = 6*m*(slope**2-residual*(1+4*qhi))
        curvature -= 4*(d-m)*(2*tolerance+arb(".18")*narrow)
        assert curvature > 17
        # If |r_i|<=eta initially, |D_i|=|r_i-mean(r)|<=4*eta/3 there.
        # The minimum has |mean(r)|<=sqrt(d/m)*eta by the trial cost.
        shift_factor = (1+ratio.sqrt())/slope
        centered_factors = [arb(4)/3+v*shift_factor for v in centered]
        mean_factor = ratio*sum(a*v for a,v in zip(centered_factors,centered))/(3*slope)
        amp = max(centered_factors)+mean_factor
        assert amp*arb("3.18e-7") < tolerance
        print("leads",count,"residual amplification below",amp.upper().str(8),
              "exterior gap above",exterior.lower().str(8))
    print("a 11431 fit at 3.18e-7 suffices for an actual gauge-matched ground")


if __name__ == "__main__":
    verify()
