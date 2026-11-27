"""第五壳完整夸克与第三枢轴五式的同一新增生成元区间证书。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


NUMBER = 2_306_888_100_090
NAME = "ch39-fifth-third-pivot-joint"


def adj(matrix):
    return matrix.conjugate().transpose()


def local_box():
    high = (".0002098522048887756", ".005627042035747022", ".005265191544896900",
            ".0002685612201948932", ".002652254968126031", ".004068660058936392",
            ".0004350425356049375", ".0006114581862599769", ".0003169594782921425",
            ".000014885208288039605", ".002373411247118261", ".0002792487339179694",
            ".08157254828963266", ".02231813110874107", ".01827618757992609", ".03551601203411357")
    return dict(layers=(9, 0, 16), centers=list(map(arb, ("4.128821611426642e-34",
                "8.737313399162484e-10", "1.668985900039446e-23"))),
                radii=list(map(arb, ("1.8e-9", "7e-9", "2e-8"))),
                high_centers=list(map(arb, high)), high_log_radius=arb(".005"))


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    ctx.prec = 384
    adjust = lib("fourth-third-pivot-three-adjusters")
    box = adjust["verify"](root, number=NUMBER, box_data=local_box())
    old = lib("fourth-third-pivot-joint")
    z0, z, dz, widths, upstream, incoming, levels = old["verify"](root, return_chart=True)
    quark = lib("fifth-two-charge-full-quark")
    shape = lib("independent-unitary-quark-right-holdout")["coordinates"]
    turn = quark["rotation"](quark["QUATERNION"])
    make_block = lib("next-coherent-full-ckm")["block_from_coordinates"]
    projectors = lib("free-four-site-native-positive-spectrum")["native"](root, NUMBER)
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    ctx.prec = 384
    z2 = z*z
    wide_ops, point_ops = (adj(z2), z), (adj(z0*z0), z0)
    varied_ops = ([adj(e*z+z*e) for e in dz], dz)
    sectors, point_sectors, tangents = [], [], []
    for iq, (q, (src, _, w, _, goal), rows) in enumerate(zip(
            (-4, 2), upstream, projectors, strict=True)):
        target = turn*make_block(goal.entries())*adj(turn)
        dst = spec(shell, q, NUMBER)[1]
        rest = (adj(dst), w, arb_mat([shape(r) for r in rows]).transpose(),
                arb_mat([[v] for v in shape(target)]))
        sectors.append((wide_ops[iq]*src, *rest))
        point_sectors.append((point_ops[iq]*src, *rest))
        tangents.append([e*src for e in varied_ops[iq]])
    inc = z2*z*incoming
    tangents.append([(e*z2+z*e*z+z2*e)*incoming for e in dz])
    source_chart = (point_sectors, z0**3*incoming), tangents, widths
    readout = adj(spec(shell, 6, NUMBER)[1])
    _, _, pivots, _ = next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root, NUMBER) if row[0] == 6)
    weights = [a.real for a in pivots[2]]
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [NUMBER])
    data = sectors, inc, readout, levels, (moments, heat, field), weights
    sites = old["solve_matrix"](root, data, box["centers"], box["radii"], name=NAME, upstream=source_chart)
    high = [x for i, x in enumerate(sites) if i not in box["layers"]]
    shifts = [abs((x/c).log()).upper() for x, c in zip(high, box["high_centers"], strict=True)]
    maximum = max(shifts)
    assert maximum < box["high_log_radius"]
    print("十六高层截止的最大对数偏移与最小截止平方", maximum.str(12),
          min(x.lower() for x in high).str(12), flush=True)
    return sites


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
