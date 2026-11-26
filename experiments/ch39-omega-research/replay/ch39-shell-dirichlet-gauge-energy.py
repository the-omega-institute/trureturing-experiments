"""Eliminate the interior of the unit-edge lifted shell with four gauge terminals."""
from pathlib import Path
from runpy import run_path

from flint import arb, fmpq, fmpq_mat

root = Path(__file__).parent
local = run_path(str(root/"ch39-local-product-contact.py"))
floating = run_path(str(root/"ch39-floating-clock-ground.py"))


def terminal_response():
    sites, pole = local["lifted_sites"]()
    laplacian, origin = local["lifted_laplacian"]()
    nearest = sorted({p for p in sites if local["geometry"]["squared"](p,pole)==68})
    boundary = [i for i,p in enumerate(sites) if i==origin or p in nearest]
    interior = [i for i in range(25) if i not in boundary]
    assert len(nearest)==3 and len(interior)==15
    def block(rows,cols):
        return fmpq_mat([[laplacian[i,j] for j in cols] for i in rows])
    inside = block(interior,interior)
    # Each leading minor is positive, so the harmonic extension uniquely minimizes.
    assert all(block(interior[:j],interior[:j]).det() > 0 for j in range(1,16))
    lift = fmpq_mat([[int(sites[i]==p) for p in nearest] for i in boundary])
    coupling = block(interior,boundary)*lift
    response = lift.transpose()*block(boundary,boundary)*lift-coupling.transpose()*inside.inv()*coupling
    common, contrast = fmpq(21,4), fmpq(315,796)
    assert all(response[i,j] == (common+2*contrast if i==j else -contrast)
               for i in range(3) for j in range(3))
    shared_excess = (laplacian[origin,origin]-3*common)/2
    assert shared_excess == fmpq(9,8)
    print("four-terminal shell: common", common, "pairwise", contrast, "shared excess", shared_excess)
    return arb(common)/2, arb(contrast)/2


def verify():
    common, pairwise = terminal_response()
    floating["verify"]()
    eps = arb("1.005e-7")
    # The extra pairwise term is nonnegative and independent of the zero label.
    extra_trial = 3*pairwise*(2*eps)**2
    assert common > arb(".5") and arb("1.55e-7")-extra_trial > arb("1.54e-7")
    ratio = (common+3*pairwise)/common
    residual_bound = (3*ratio).sqrt()*eps
    radius = (2*residual_bound+2*eps)/arb(".1435")
    target = floating["gauge"]["shortest"]["model"]["target"]
    center = target(11431,arb("26.81676714125598"))
    low,high = center*(-2*radius).exp(), center*(2*radius).exp()
    slope = 2*low-arb(".09")
    residual = eps+arb(".09")*radius+high-low
    curvature = 6*common*(slope**2-residual*(1+4*high))
    # Each difference of two gauge accelerations is bounded in modulus by two.
    spread = 2*eps+arb(".18")*radius
    curvature -= 12*pairwise*spread
    # At the minimum the mean residual is <= sqrt(ratio)*eps and increases
    # at least at the displayed slope throughout the possible minimum interval.
    location = (1+ratio.sqrt())*eps/slope
    assert curvature > 660 and location < arb("3.26e-8")
    assert residual_bound < arb("1.94e-7")
    # The weights below are independent of scale because every target has
    # the same exponential factor. They rule out affine equilibrium sources.
    time = arb("26.81676714125598")
    lower, upper = target(2994,time), target(14840,time)
    weight = (center-upper)/(lower-upper)
    assert 0 < weight < 1
    free = floating["penalty"]["doublet"]["joint"]["contact"]["clock"]["free_ground"]
    score = floating["penalty"]["competition"]["score"]
    gain = free(score(11431),2)-weight*free(score(2994),2)-(1-weight)*free(score(14840),2)
    assert gain > arb(".0000088")
    print("Dirichlet minimum scale radius",location.upper().str(9),
          "curvature above",curvature.lower().str(9),"max residual",residual_bound.upper().str(9))


if __name__ == "__main__":
    verify()
