"""用同一个四站嵌入，在中央和右邻实际壳分别重建六个汤川中心与完整 CKM。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


CENTRAL = 2_306_888_100_085


def adj(matrix):
    return matrix.conjugate().transpose()


def common_embedding(u_text=".00384", v_text=".0384", quaternion=(".15495", ".07140", ".97665", "-.13053")):
    u, v = arb(u_text), arb(v_text)
    a, b = (u*v).sqrt(), (1-u*v).sqrt()
    p = v/(u+v)
    f, g = p.sqrt(), (1-p).sqrt()
    intrinsic = acb_mat([[1, 0, 0], [0, 0, b], [0, g, -a*f], [0, -f, -a*g]])
    seed = [arb(s) for s in quaternion]
    length = sum(s*s for s in seed).sqrt()
    x, y, z, w = [s/length for s in seed]
    h = 1/arb(2).sqrt()
    rotation = acb_mat([[acb(w, -z), acb(-y, -x), 0, 0],
                        [acb(y, -x), acb(w, z), 0, 0],
                        [0, 0, h, h], [0, 0, acb(0, h), acb(0, -h)]])
    return u, v, intrinsic, rotation


def extend_target(matrix, u, v):
    assert matrix[0, 1].contains(0) and matrix[1, 2].imag.contains(0)
    d, s, t = [matrix[j, j].real for j in range(3)]
    k, w = matrix[0, 2], matrix[1, 2].real
    b = (1-u*v).sqrt()
    low = (t-u*v*s-(u-v)*w)/(1-u*v)
    plus, minus = s-w/u, s+w/v
    return acb_mat([[d, k/b, 0, 0], [k.conjugate()/b, low, 0, 0],
                    [0, 0, plus, 0], [0, 0, 0, minus]])


def verify_nodes(root, nodes, basis_coefficient, geometry, chosen):
    single = run_path(str(root/"ch39-complex-four-site-unit-quark-amplitudes.py"))
    fixed = run_path(str(root/"ch39-free-four-site-native-positive-spectrum.py"))
    ckm = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    flavor = run_path(str(root/"ch39-three-gram-source-flavor.py"))
    diag = single["diag"]
    ctx.prec = 320
    angles = [arb(i).sin() for i in (".22704", ".04275", ".00377")]
    mixing = acb_mat(ckm["standard_matrix"](*angles, arb("1.139")))
    masses = [list(map(arb, row)) for row in ((".00000615", ".00311", ".8616"),
                                              (".0000135", ".000268", ".01401"))]
    up, down = [diag([y/row[-1] for y in row]) for row in masses]
    up = adj(mixing)*up*mixing
    basis = fixed["target_basis"](up, down*arb(basis_coefficient))*diag([acb(0, 1), 1, acb(0, -1)])
    assert (adj(basis)*basis).contains(diag([1]*3))
    original = [adj(basis)*m*basis*row[-1] for m, row in zip((up, down), masses, strict=True)]
    u, v, intrinsic, rotation = common_embedding(*geometry)
    embedding = rotation*intrinsic
    assert (adj(rotation)*rotation).contains(diag([1]*4))
    assert (adj(embedding)*embedding).contains(diag([1]*3))
    expected_j = (mixing[0, 0]*mixing[1, 1]*mixing[0, 1].conjugate()*mixing[1, 0].conjugate()).imag
    central_effective = []
    responses = {}
    for node, phase_coordinate in nodes:
        phase = acb(1, arb(phase_coordinate))
        transform = diag([phase/phase.conjugate(), 1, 1])
        assert (adj(transform)*transform).contains(diag([1]*3))
        sectors = fixed["native"](root, node)
        eigenspaces = []
        for qi, charge in enumerate((-4, 2)):
            rows = sectors[qi]
            target = transform*original[qi]*adj(transform)
            four = extend_target(target, u, v)
            assert (adj(intrinsic)*four*intrinsic-target).contains(acb_mat(3, 3))
            picks = [j-1 for j in chosen[qi]]
            _, weights = single["reconstruct"](rows, four, rotation, picks)
            responses[node,charge] = weights
            # 在原生四站坐标中叠加十九个实际谱投影，再用两壳共享的物理嵌入压缩。
            raw = sum((w*p for w, p in zip(weights, rows, strict=True)), acb_mat(4, 4))
            effective = adj(embedding)*raw*embedding
            assert (effective-target).contains(acb_mat(3, 3))
            assert max(abs(z).upper() for z in (effective-target).entries()) < arb("1e-40")
            levels, projectors = flavor["projections"](effective*effective)
            assert all(l.overlaps(y*y) for l, y in zip(levels, masses[qi], strict=True))
            eigenspaces.append(projectors)
            largest = max(w.lower() for w in weights).union(max(w.upper() for w in weights))
            print(node, charge, "四个自由层", *[weights[i].str(13) for i in picks],
                  "两个字符", four[2, 2].real.str(13), four[3, 3].real.str(13),
                  "最大", largest.str(15))
            if node == CENTRAL:
                central_effective.append(effective)
            elif central_effective:
                difference = abs(effective[0, 2]-central_effective[qi][0, 2])
                assert difference > 0
                print("两壳有效矩阵的 13 元之差的模", charge, difference.str(15))
        pu, pd = eigenspaces
        assert all((pu[i]*pd[j]).trace().real.overlaps(abs(mixing[i, j])**2)
                   for i in range(3) for j in range(3))
        jarlskog = (pu[0]*pd[0]*pu[1]*pd[1]).trace().imag
        assert jarlskog > 0 and jarlskog.overlaps(expected_j)
        print(node, "有符号 J", jarlskog.str(15))
    return embedding, responses, rotation


def verify(root):
    embedding, responses, rotation = verify_nodes(root, ((CENTRAL, "0"), (CENTRAL+1, ".0294")),
                                                 ".1", (), ((6, 7, 10, 13), (6, 7, 10, 12)))

    cover = run_path(str(root/"ch39-neighbour-all-frame-support.py"))
    cover["verify"](root, -1, rotation)
    left_geometry = (".0182", ".113", ("-.21668", "-.53516", ".71837", "-.38807"))
    _, _, left_rotation = verify_nodes(root, ((CENTRAL-1, "0"),), ".74", left_geometry,
                                       ((1, 10, 16, 18), (6, 7, 10, 18)))
    for offset in (0, 1):
        cover["verify"](root, offset, left_rotation, left_geometry[:2])
    return embedding,responses


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
