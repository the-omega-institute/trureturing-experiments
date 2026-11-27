"""十九个原生站位中混合十三站，在右外壳实现三个上型汤川耦合。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


SITES = (0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 13, 14, 17)
SCALE = 10_000
# 反厄米矩阵 K 的上三角，从对角元起依次给出高斯整数的实部和虚部。
UPPER = (
    ((0, -2850), (-1642, 1140), (-1240, 1751), (2493, -973), (-379, 247), (-912, 1097), (-70, -236), (33, -51), (-179, -54), (1146, -605), (-234, -159), (-24, -24), (842, 813),),
    ((0, -346), (-1177, -479), (-447, 549), (1109, -965), (-704, 607), (-27, -144), (-84, -168), (-34, -205), (-103, 403), (245, -39), (8, 14), (-311, -120),),
    ((0, 352), (-2445, 1045), (8, -1793), (1308, 956), (13, -142), (-87, -177), (-7, -247), (-814, -106), (-37, 85), (-2, 64), (-260, -684),),
    ((0, -123), (-1982, 1611), (1876, -933), (-163, -146), (22, -180), (-41, -3), (-358, -131), (183, 270), (18, 54), (259, -399),),
    ((0, 2317), (-1277, -2332), (-39, -111), (-85, -106), (-66, -148), (805, 935), (346, -305), (64, -22), (-499, 712),),
    ((0, 2524), (-8, -108), (-78, -175), (15, -202), (-18, -468), (-3, -222), (-40, 30), (202, -357),),
    ((0, -373), (-12, -386), (5, -390), (-156, 223), (-216, 299), (-206, -67), (-120, 313),),
    ((0, -365), (-12, -386), (-225, 306), (-90, 314), (-205, -65), (-150, 193),),
    ((0, -374), (-119, 325), (-140, 197), (-204, -68), (-240, 279),),
    ((0, -283), (1, -95), (118, -3), (163, -189),),
    ((0, -18), (123, -9), (66, -57),),
    ((0, -12), (-126, -2),),
    ((0, -207),),
)


def local_transport(adj, eye):
    n = len(SITES)
    k = acb_mat(n, n)
    for i, row in enumerate(UPPER):
        assert len(row) == n - i and row[0][0] == 0
        for j, (real, imag) in enumerate(row, start=i):
            k[i, j] = acb(real, imag)
            k[j, i] = acb(-real, imag)
    assert k == -adj(k)
    scaled_identity = eye(n) * SCALE
    block = (scaled_identity + k).solve(scaled_identity - k)
    assert (adj(block) * block).contains(eye(n))
    assert all(any(abs(block[i, j]) > 0 for j in range(n) if i != j) for i in range(n))
    transport = eye(19)
    for i, physical_i in enumerate(SITES):
        for j, physical_j in enumerate(SITES):
            transport[physical_i, physical_j] = block[i, j]
    return transport


def verify(root):
    training = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    transport = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    up = run_path(str(root / "ch39-rank-two-generator-up-quark.py"))
    source_node = training["CENTRAL"] + 1
    _, response, _ = training["verify_nodes"](
        root, ((source_node, ".0294"),), ".1", (), ((6, 7, 10, 13), (6, 7, 10, 12))
    )
    ctx.prec = 384
    adj = up["adj"]
    unitary = local_transport(adj, up["eye"])
    _, source = transport["spectrum"](shell, -4, source_node)
    _, target = transport["spectrum"](shell, -4, source_node + 1)
    weights = up["read"](source, target, response[source_node, -4], unitary)
    sites = (shell["soft"]["pole"], *shell["boundary_sites"]())
    assert sites == (14, 8, 7, 6)
    c = acb_mat([[target[i, j] for j in range(19)] for i in sites])
    weighted = acb_mat([[c[i, j] * weights[j] for j in range(19)] for i in range(4)])
    print("活动站位", SITES)
    return up["compression"](weighted * adj(c), list(map(arb, up["CENTERS"])))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
