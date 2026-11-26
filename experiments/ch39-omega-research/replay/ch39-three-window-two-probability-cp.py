"""三窗零点同时固定两个夸克跃迁概率时，夹定平均 CP 亏额的两端。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

previous = run_path(str(Path(__file__).with_name("ch39-three-window-extreme-cp.py")))
cp = previous["cp"]
outer_interval = (".97561882", ".998810426891862")
inner_centers = [
    (".15674163", "-.05916834"), (".072", ".072"),
    ("-.01585316297", "-.12259306920"), (".12240241897", ".01575652959"),
    (".03410521449", "-.15106657023"), (".03408818760", "-.15106184447"),
]


def observables():
    ctx.prec, ctx.cap = 384, 3
    flavor = cp["flavor"]
    data = cp["two"]["spectral_data"]()
    polys = cp["cp_data"](data, flavor["reference_matrix"]())
    ref = flavor["singlet"](11661)[1]
    up_ref, charm_ref = cp["two"]["probabilities"](ref, data)
    current_ref = cp["cp_value"](ref, polys)
    slope = charm_ref[1]/up_ref[1]

    def jets(c):
        up, charm = cp["two"]["probabilities"](c, data)
        excess = (charm-charm_ref[0]-slope*(up-up_ref[0]))*10**12
        loss = 1-cp["cp_value"](c, polys, jet=True)/current_ref
        return (up-up_ref[0])*10**5, excess, loss

    return jets


def verify_plane(jets, coefficients, direction):
    intercept, up_slope, charm_slope = [arb(x) for x in coefficients]

    def difference(values, order):
        up, charm, loss = values
        value = loss[order] - up_slope*up[order] - charm_slope*charm[order]
        return direction*(value-intercept if order == 0 else value)

    pending = [tuple(map(arb, outer_interval))]
    while pending:
        left, right = pending.pop()
        mid = (left+right)/2
        width = right-mid
        point = jets(mid)
        region = jets(left.union(right))
        residual = (difference(point, 0) + difference(point, 1)*arb(0, width)
                    + difference(region, 2)*arb(0).union(width*width))
        if residual > 0:
            continue
        assert width > arb("1e-10")
        pending.extend(((left, mid), (mid, right)))
    return intercept + 2*charm_slope


def box_values(jets, c):
    mid = c.mid()
    return [point[0] + box[1]*(c-mid) for point, box in zip(jets(mid), jets(c), strict=True)]


def pair_at_reference(first, second):
    weight = -second[0]/(first[0]-second[0])
    assert 0 < weight < 1
    return [weight*x + (1-weight)*y for x, y in zip(first[1:], second[1:], strict=True)]


def match_both(endpoint, left, right):
    charm_a, loss_a = pair_at_reference(endpoint, left)
    charm_b, loss_b = pair_at_reference(endpoint, right)
    weight = (2-charm_b)/(charm_a-charm_b)
    assert 0 < weight < 1
    return weight*loss_a + (1-weight)*loss_b


def verify():
    coords = previous["open_windows"](inner_centers)
    jets = observables()
    flavor = cp["flavor"]
    axes = [flavor["singlet_at"]({724: z, 804: acb(1)})[1] for z in (acb(1), acb(0, 1))]
    lo, hi = map(arb, outer_interval)
    assert axes[0]*arb(previous["aperture"]).cos()-abs(axes[1]) > lo
    assert sum(v*v for v in axes) < hi*hi
    lower = verify_plane(jets, ("-.01495012", "-.47587894", ".03198095"), 1)
    upper = verify_plane(jets, (".02143693", "-.51013474", ".03885414"), -1)
    points = [box_values(jets, c) for c in coords]
    near_low = match_both(points[0], points[2], points[3])
    near_high = match_both(points[1], points[4], points[5])
    assert lower < near_low < arb(".04901458")
    assert arb(".09914505") < near_high < upper
    print("两个精确概率下的外下界、可达低值、可达高值、外上界",
          *[v.str(13) for v in (lower, near_low, near_high, upper)])


if __name__ == "__main__":
    verify()
