"""在第五壳实际四站上检验旧四阶目标，认证两荷独立输运的共同三代压缩。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


NUMBER = 2_306_888_100_090
QUATERNION = (30903302781, 11011638290, 126118559266, 991472560744)
SEEDS = (
    (None, 947553943479, 947553943479, None, 3382481944, 1919765509, 14154325671,
     936656803820, 1919765509, 600363689, None, 12910925441, 1919765509,
     947553943479, 149010239048, None, None, None, None),
    (None, 219791080379, 46954358380, 153774785360, 153774785360, 74636520795, None,
     153774785360, 113277211006, 10452567, None, 4092715845, 79761243806,
     153774785360, 10453497, None, None, None, None),
)


def rotation(quaternion):
    length = arb(sum(v*v for v in quaternion)).sqrt()
    x, y, z, w = [arb(v)/length for v in quaternion]
    return acb_mat([[acb(w, -z), acb(-y, -x), 0, 0], [acb(y, -x), acb(w, z), 0, 0],
                    [0, 0, 1, 0], [0, 0, 0, 1]])


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    independent = lib("independent-unitary-quark-right-holdout")
    emb = lib("common-embedding-central-right-ckm")
    adj = emb["adj"]
    coords = independent["coordinates"]
    original = lib("one-generator-gravity-full-ckm")
    reference = original["targets"](root, independent, emb)
    source = independent["CENTER"]+1
    _, training, _ = emb["verify_nodes"](root, ((source, ".0294"),), ".1", (),
                                         ((6, 7, 10, 13), (6, 7, 10, 12)))
    sectors = lib("free-four-site-native-positive-spectrum")["native"](root, NUMBER)
    ctx.prec = 384
    old_rotation = rotation(lib("fourth-full-ckm-gravity")["QUATERNION"])
    turn = rotation(QUATERNION)
    _, _, intrinsic, initial = emb["common_embedding"](*independent["GEOMETRY"])
    frame = turn*old_rotation*initial*intrinsic
    eye = lib("rank-two-generator-up-quark")["eye"]
    assert (adj(frame)*frame-eye(3)).contains(acb_mat(3, 3))
    make_block = lib("next-coherent-full-ckm")["block_from_coordinates"]
    margins = []
    for i, (charge, rows) in enumerate(zip((-4, 2), sectors, strict=True)):
        known = training[source, charge]
        assert all(v > 0 for v in known)
        ref = make_block(reference[i].entries())
        old = old_rotation*ref*adj(old_rotation)
        if i == 0:
            # 重排不等式给出这个实线性泛函在整个 Schur–Horn 多面体上的上界。
            separator = (-11010075899, -1734690587, -11545700384, -9223656529, 0, 0)

            def value(matrix):
                return sum(q*x for q, x in zip(separator, coords(matrix), strict=True))/100000000

            # 各阶训练振幅以独立排序的上下端点包围，包括重特征值。
            bounded = [a.union(b) for a, b in zip(sorted(v.lower() for v in known),
                                                  sorted(v.upper() for v in known), strict=True)]
            upper = sum((a*b).upper() for a, b in zip(sorted(value(row).upper() for row in rows),
                                                      bounded, strict=True))
            gap = value(old)-upper
            print("旧上型四阶目标与可达集的分离量下界", gap.lower().str(15), flush=True)
            assert gap > arb(".0052")
        goal = turn*old*adj(turn)
        block, margin = independent["reconstruct"](rows, goal, known, SEEDS[i])
        diff = adj(frame)*block*frame-adj(initial*intrinsic)*ref*initial*intrinsic
        error = max(abs(v).upper() for v in diff.entries())
        assert diff.contains(acb_mat(3, 3)) and error < arb("1e-38")
        assert margin > arb(("6.58e-8", "1.46e-8")[i])
        print("荷、十八个主序裕量的最小下界、三代矩阵区间误差",
              charge, margin.str(14), error.str(8), flush=True)
        margins.append(margin)
    return margins


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
