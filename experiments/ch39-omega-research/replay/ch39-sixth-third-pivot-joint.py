"""第六壳的一个新增酉同时实现两荷完整三代目标和第三枢轴五个端点条件。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


NUMBER = 2_306_888_100_091
NAME = "ch39-sixth-third-pivot-joint"
QUATERNION = (-2465348963, -23139382617, -109355319385, 993730297993)


def adj(matrix):
    return matrix.conjugate().transpose()


def local_box():
    high = ("4.9635728605016545e-05", ".006273050296963615", ".004077759097579717",
            ".0012877458874353703", ".0017651988249919688", ".0007514634181214349",
            ".0036434314175488054", ".0006751158881955718", "3.868194113568315e-05",
            "4.7668953895416965e-05", ".005166719592797377", ".000964916575092361",
            ".045583322428706864", ".024845217554012533", ".01993238726680459", ".046192578029678764")
    return dict(layers=(9, 0, 16), centers=list(map(arb, ("3.914491233849589e-34",
                "7.480893435045741e-10", "8.849433273594449e-23"))),
                radii=list(map(arb, ("1.8e-9", "1.0e-8", "2e-8"))),
                high_centers=list(map(arb, high)), high_log_radius=arb(".003"))


def inputs(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    generator, ancestor = lib("fifth-third-pivot-joint")["verify"](root, return_chart=True)
    w0, w, dw, widths, previous, incoming, levels = generator
    (point_previous, point_inc), dz, z_widths = ancestor
    qlib = lib("fifth-two-charge-full-quark")
    shape = lib("independent-unitary-quark-right-holdout")["coordinates"]
    make_block = lib("next-coherent-full-ckm")["block_from_coordinates"]
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    projectors = lib("free-four-site-native-positive-spectrum")["native"](root, NUMBER)
    ctx.prec = 384
    turn = qlib["rotation"](QUATERNION)
    w2 = w*w
    ops, refs = (adj(w2), w), (adj(w0*w0), w0)
    varied = ([adj(d*w+w*d) for d in dw], dw)
    sectors, point_sectors, tangents = [], [], []
    for iq, (q, (src, _, weights, _, goal), rows) in enumerate(zip(
            (-4, 2), previous, projectors, strict=True)):
        target = turn*make_block(goal.entries())*adj(turn)
        dst = spec(shell, q, NUMBER)[1]
        rest = (adj(dst), weights, arb_mat([shape(r) for r in rows]).transpose(),
                arb_mat([[v] for v in shape(target)]))
        sectors.append((ops[iq]*src, *rest))
        point_sectors.append((refs[iq]*point_previous[iq][0], *rest))
        tangents.append([d*src for d in varied[iq]]+[ops[iq]*d for d in dz[iq]])
    w3 = w2*w
    inc = w3*incoming
    tangents.append([(d*w2+w*d*w+w2*d)*incoming for d in dw]+[w3*d for d in dz[2]])
    chart = (point_sectors, w0**3*point_inc), tangents, widths+z_widths
    readout = adj(spec(shell, 6, NUMBER)[1])
    _, _, pivots, _ = next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root, NUMBER) if row[0] == 6)
    weights = [a.real for a in pivots[2]]
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [NUMBER])
    return (sectors, inc, readout, levels, (moments, heat, field), weights), chart


def solve(root, box, data, chart):
    matrix = run_path(str(root/"ch39-fourth-third-pivot-joint.py"))["solve_matrix"]
    sites = matrix(root, data, box["centers"], box["radii"], name=NAME, upstream=chart, component_box=True)
    high = [x for i, x in enumerate(sites) if i not in box["layers"]]
    maximum = max(abs((x/c).log()).upper() for x, c in zip(high, box["high_centers"], strict=True))
    assert maximum < arb("6.1e-6") < box["high_log_radius"]
    assert min(x.lower() for x in high) > arb("3.8681e-5")
    print("十六高层截止的最大对数偏移与最小截止平方", maximum.str(12),
          min(x.lower() for x in high).str(12), flush=True)
    return sites


def verify(root):
    ctx.prec = 384
    adjust = run_path(str(root/"ch39-fourth-third-pivot-three-adjusters.py"))
    box = adjust["verify"](root, number=NUMBER, box_data=local_box())
    return solve(root, box, *inputs(root))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
