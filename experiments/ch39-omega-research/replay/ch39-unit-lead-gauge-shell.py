"""Couple four external source voltages by unit edges to a fully free lifted shell."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, fmpq, fmpq_mat

previous = run_path(str(Path(__file__).with_name("ch39-shell-dirichlet-gauge-energy.py")))
local, floating = previous["local"], previous["floating"]


def terminal_modes(count):
    sites,pole = local["lifted_sites"]()
    matrix,origin = local["lifted_laplacian"]()
    nearest = sorted({p for p in sites if local["geometry"]["squared"](p,pole) == 68})
    leads = [[origin]]+[[i for i,v in enumerate(sites) if v==p][:count] for p in nearest]
    inside, coupling, outer = fmpq_mat(matrix), fmpq_mat(25,4), fmpq_mat(4,4)
    for j,group in enumerate(leads):
        for i in group:
            inside[i,i] += 1
            coupling[i,j] -= 1
            outer[j,j] += 1
    response = outer-coupling.transpose()*inside.inv()*coupling
    assert all(sum(response[i,j] for j in range(4))==0 for i in range(4))
    common = (response[1,1]+2*response[1,2])/2
    contrast = (response[1,1]-response[1,2])/2
    assert all(response[i,j]/2 == (common-contrast)/3+(contrast if i==j else 0)
               for i in range(1,4) for j in range(1,4))
    return common, contrast


def verify(*, target_error=0):
    floating["verify"]()
    expected = [(fmpq(7,61),fmpq(427,1129)),(fmpq(28,209),fmpq(1708,2381)),
                (fmpq(21,148),fmpq(1281,1252))]
    for count,case in enumerate(expected,1):
        assert terminal_modes(count)==case
    # The reverse direction already has a global ceiling. Evaluate the
    # remaining forward scores at the additional indexed zeros on both sides.
    pulse = floating["ceiling"]["pulse"]
    fringe_ceiling = arb(".00741904")
    for first,stop in ((10400,10815),(12093,12500)):
        for start in range(first,stop,32):
            for z in acb.zeta_zeros(start,min(32,stop-start)):
                a,b = pulse["port_outputs"](pulse["generator"](z.imag),2,impulse=True)
                assert abs((-7*a[0]+24*b[0])/25)**2 < fringe_ceiling
    pen = floating["penalty"]
    ref_score = pen["competition"]["score"](11431)
    assert ref_score > arb(".00742343")
    # A mixed orientation has no greater score than its better constituent.
    assert arb(".007300684") < fringe_ceiling
    free = pen["doublet"]["joint"]["contact"]["clock"]["free_ground"]
    fringe_gap = free(fringe_ceiling,2)-free(ref_score,2)
    assert fringe_gap > arb("8e-6")

    t = arb("26.81676714125598")
    # Allow monotone substitutes with the same exponential scale factor.
    error = arb(target_error)
    assert error >= 0
    eps = arb("1.005e-7")+error
    target = floating["gauge"]["shortest"]["model"]["target"]
    center = target(11431,t)+arb(0,error.upper())
    radius = arb(".061")
    edge = min((target(10399,t+radius)-center).lower(),(center-target(12500,t-radius)).lower())
    edge -= error+eps+arb(".09")*radius
    spread = arb(".1435")*radius-2*eps
    assert edge>0 and spread>0
    common = arb(min(m for m,d in expected))
    contrast = arb(min(d for m,d in expected))
    maximum = arb(max(d for m,d in expected))
    trial = maximum*3*eps**2
    exterior = min((contrast*spread**2/2).lower(),(3*common*edge**2).lower())-trial-arb(".000011193")
    assert exterior > arb("3e-6") and fringe_gap-trial > arb("8e-6")
    ratio = maximum/common
    narrow = (2*(3*ratio).sqrt()+2)*eps/arb(".1435")
    low,high = center*(-2*narrow).exp(),center*(2*narrow).exp()
    slope = 2*low-arb(".09")
    r = eps+arb(".09")*narrow+high-low
    curvature = 6*common*(slope**2-r*(1+4*high))-4*(maximum-common)*(2*eps+arb(".18")*narrow)
    location = (1+ratio.sqrt())*eps/slope
    assert curvature>28 and location < arb("6.2e-8") and (3*ratio).sqrt()*eps < arb("5.3e-7")
    # Sharpen the readout at the stationary scale in each individual network.
    for m,d in expected:
        # Differentiate 3*m*mean_residual^2 + d*sum_i(X_i-Xbar)^2.
        # Centered gauge values at t are bounded by 2eps; centered speeds by .18.
        mean_bound = arb(d)/arb(m)*(2*eps+arb(".18")*location)*arb(".18")/slope
        max_error = mean_bound + 2*eps+arb(".18")*location
        assert max_error < arb("2.6e-7") < arb(2)**-21
        print("unit-lead modes",m,d,"three residuals below",max_error.upper().str(7))
    print("one to three leads per gauge channel: margin above",exterior.lower().str(8),
          "scale radius",location.upper().str(8))


if __name__ == "__main__":
    verify()
