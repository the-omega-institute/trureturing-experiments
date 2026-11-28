"""验证两枚实际零点的共同可加谱作用需要至少十二次。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


# 相位均已乘部门荷。前一组相对于中央，后一组是绝对相位，单位为弧度。
PHASES = {
    -4: (
        (("0", "0"), ("0", "-1e-8"), ("1e-6", "1e-6")),
        (("-.61536", "-.16194"), ("-.6852883", "-1.1298121"), (".6253561", "1.1579139")),
    ),
    2: (
        (("0", "0"), ("1e-8", "1e-8"), ("-1e-6", "0"), ("-1e-6", "-1e-6")),
        (("-1.06879", "-1.13356"), ("-.2266231", "-1.4679559")),
    ),
}


def lowest(matrix):
    spectrum = matrix.eig(algorithm="rump")
    assert all(z.is_finite() and z.imag.contains(0) for z in spectrum)
    levels = sorted((z.real for z in spectrum), key=lambda x: float(x.mid()))
    assert all(levels[0] < x for x in levels[1:])
    return levels[0]


def verify(replay):
    selector = run_path(str(replay / "ch39-light-spectrum-zero-selector.py"))
    shell = run_path(str(replay / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay / "ch39-native-light-projector-pump.py"))
    ctx.prec = 384
    blocks = list(pump["rotation_bases"](shell))
    adj = pump["adjoint"]
    moments = selector["chebyshev_moments"]
    vectors = []
    for charge, index in ((-4, 2345666), (2, 964890)):
        height = acb.zeta_zero(index).imag
        central_matrix = shell["laplacian"](charge, height)
        central = moments([lowest(adj(b)*central_matrix*b) for b in blocks], 11)
        phases = shell["geo"]["phases"](height)
        center_factors = [phases[d]**charge for d in (724, 804)]
        base, edges = pump["clock_matrices"](shell, charge)
        relative, absolute = PHASES[charge]
        for origin, offsets in ((center_factors, relative), ((1, 1), absolute)):
            for point in offsets:
                factors = [z*acb(0, arb(t)).exp() for z, t in zip(origin, point, strict=True)]
                link = edges[0]*factors[0] + edges[1]*factors[1]
                flat = base + link + adj(link)
                candidate = moments([lowest(adj(b)*flat*b) for b in blocks], 11)
                displacement = [v-w for v, w in zip(candidate, central, strict=True)]
                length = sum(v**2 for v in displacement).sqrt()
                assert length > 0
                vectors.append([v/length for v in displacement])
    radius = selector["interior_radius"](vectors)
    assert radius > arb("2.75e-12")
    print("十二维增广矩阵给出的十一维归一位移内半径", radius.str(16))


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
