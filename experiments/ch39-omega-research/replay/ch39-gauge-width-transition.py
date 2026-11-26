"""Bound the input-box widening that changes the selected arithmetic height."""
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx

root = Path(__file__).parent
model = run_path(str(root/"ch39-floating-scale-winner.py"))
hardy = run_path(str(root/"ch39-hardy-height-selection.py"))
q = model["q"]


def robust_exclusion(delta=q(".001007"), time=q("26.861145"), cutoff=10603):
    full = (1-delta).union(1+delta)
    corners = [1-delta,1-delta+q(".0002"),1+delta-q(".0002"),1+delta]
    low,weak_cut,color_cut,high = corners
    faces = [(weak_cut.union(high),full),(low.union(weak_cut),low.union(color_cut)),
             (low.union(weak_cut),color_cut.union(high))]
    offset = q(model["fixed"]["OFFSET"])
    reference = acb.zeta_zero(10854).imag
    target = model["target"](10854,time)*(reference**2+offset)/(q(cutoff)**2+offset)
    logs = [2*(v*q(".99").union(q("1.01"))).log() for v in model["yukawa_centers"]]
    margins = []
    for weak,color in faces:
        inverse = [1/(v*f)**2 for v,f in zip(model["centers"],[full,weak,color])]
        end,_,diff = model["evolve"](inverse,logs,time,16384)
        margins.append((target-end[1]/model["ratios"][1],-diff))
    assert min(x.lower() for x,_ in margins)>q(".00004")
    assert min(x.lower() for _,x in margins)>q(".000006")
    def monotonic(step, inverse, rate):
        if step >= 23*64:
            weak = rate[1]/model["ratios"][1]
            assert 0 < weak < rate[2]
    model["running"]["integrate"]([v*full for v in model["centers"]],33*64,monotonic)
    print("Gauge exclusion up to",cutoff,"at input half-width",delta.str(7),"margins",
          min(x.lower() for x,_ in margins).str(8),min(x.lower() for _,x in margins).str(8))


def widened_witness():
    time = q("26.861182485")
    common = model["target"](10854,time)
    logs = list(map(q, ("-25.4446915838 -12.9928374813 -1.5528331663 -23.8189018241 "
                        "-17.8423094899 -10.1208976879 -25.6505068564 -14.9499484326 "
                        "-9.2850641599").split()))
    inv,y,_ = model["evolve"]([common*r for r in model["ratios"]],logs,-time,131072,
                               diagonal_quarks=True)
    shifts = [(1/(x.sqrt()*c)-1)*10**6 for x,c in zip(inv,model["centers"])]
    assert max(abs(x).upper() for x in shifts)<1011
    assert all(q(".99") < (v/2).exp()/c < q("1.01") for v,c in zip(y,model["yukawa_centers"]))
    earlier,old = [acb.zeta_zero(n).imag for n in (10854,10875)]
    gain = hardy["weight"](earlier)/hardy["weight"](old)-1
    assert gain > q(".0001")
    # A nearby 64-bit grid point also fits: the same endpoint remains exact after retiming.
    grid = 2**64
    offset = q(model["fixed"]["OFFSET"])
    finite = run_path(str(root/"ch39-finite-arithmetic-height.py"))
    def grid_score(height):
        return hardy["weight"](height)/(1+abs(finite["finite_value"](height))/256)
    old_grid = arb(195923358475422178854846)/grid
    for numerator in [(earlier*grid).floor(),(earlier*grid).ceil(),(earlier*grid).ceil()+1]:
        x = numerator/grid
        retiming = ((earlier**2+offset)/(x**2+offset)).log()/2
        assert abs(retiming)<q("1e-22")
        grid_inv,grid_y,_ = model["evolve"](inv,y,-retiming,1,diagonal_quarks=True)
        assert all(abs(1/(u.sqrt()*c)-1)<q(".001011") for u,c in zip(grid_inv,model["centers"]))
        assert all(q(".99") < (v/2).exp()/c < q("1.01") for v,c in zip(grid_y,model["yukawa_centers"]))
        assert grid_score(x)>grid_score(old_grid)*(1+q(".00160"))
    print("1011 ppm witness for zero 10854, initial gauge shifts",*[a.str(8) for a in shifts],
          "relative gain",gain.str(9))


if __name__ == "__main__":
    ctx.prec = 144
    robust_exclusion()
    widened_witness()
