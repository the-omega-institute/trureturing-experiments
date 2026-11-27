"""第五壳给定两荷输出列时，六维补空间的共同第四步能否线性唯一恢复。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


PROBES = ([(-13859, 88646), (-22272, 82682), (-11813, 101829), (-21759, 75356), (-21324, 67151), (-5342, 70577),
  (-55073, 220408), (-57023, 233084), (-100601, 213095), (35877, 20683), (33159, 21076),
  (17766, 40887), (9904, 47462), (48431, 12093), (-849515, -238546), (28779, 7091), (36300, 12373),
  (29956, 12215), (2669, 34477)],
 [(-84249, -102549), (-38927, -11648), (-156914, -64366), (-133904, 12277), (-22175, -45569),
  (-77077, -45346), (-43276, 76554), (-58173, -40095), (-102514, 11730), (-119918, -152246),
  (12307, 79533), (-158447, -251238), (-82737, 148751), (-263159, 91689), (-277650, -525853),
  (-97981, 261398), (154957, -135231), (137820, -402571), (74195, 25627)])


def adj(matrix):
    return matrix.conjugate().transpose()


def positive(matrix):
    size = matrix.nrows()
    lower = acb_mat(size, size)
    for i in range(size):
        pivot = matrix[i, i] - sum(abs(lower[i, k])**2 for k in range(i))
        assert pivot.imag.contains(0) and pivot.real > 0
        lower[i, i] = pivot.real.sqrt()
        for j in range(i + 1, size):
            value = matrix[j, i] - sum(lower[j, k]*lower[i, k].conjugate() for k in range(i))
            lower[j, i] = value / lower[i, i]


def lower_mass_outside(step, source, levels, sectors, high, probes):
    low, heavy, _ = levels
    hminus = heavy - low
    for index, (pairs, bound) in enumerate(zip(probes, (".206", ".042"), strict=True), 1):
        src, _, weights, _, _ = sectors[1]
        omitted = [i for i, w in enumerate(weights) if w.overlaps(weights[index])]
        e = acb_mat([[src[i,j] for j in omitted] for i in range(19)])
        q = acb_mat([[acb(a,b)] for a,b in pairs])
        q = q / arb(sum(a*a+b*b for a,b in pairs)).sqrt()
        def upper_norm(row):
            return sum(abs(a).upper()**2 for a in row.entries()).sqrt().upper()
        old_row = adj(q) * step**3 * source
        observed = sum(max(arb(0), abs(a).lower())**2 for a in old_row.entries()).sqrt().lower()
        leakage = 2*sum(upper_norm(adj(q)*step**j*e) for j in range(1,4))
        head = upper_norm(adj(q)*high)
        gap = observed - leakage - head
        assert gap > 0, (observed, leakage, head)
        residual = hminus*gap**2
        print("自由六重簇代表层、旧读数、三步误差、高三层投影、外部振幅下界",
              index+1, *[v.str(10) for v in (observed, leakage, head, residual)], flush=True)
        assert residual > arb(bound)


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    center, error, sectors, source, levels = lib("fifth-third-pivot-joint")["verify"](
        root, return_transport=True
    )
    size = 19
    radius = arb(0, error.upper())
    step = acb_mat([[center[i, j] + acb(radius, radius) for j in range(size)] for i in range(size)])
    eye = lib("rank-two-generator-up-quark")["eye"]
    frames = []
    for src, _, weights, _, _ in sectors:
        # 源训练把同一字符的六层赋同一标量；区间分离另外十三层。
        omitted = [i for i, value in enumerate(weights) if value.overlaps(weights[2])]
        kept = [i for i in range(size) if i not in omitted]
        assert len(omitted) == 6 and len(kept) == 13
        def columns(indices):
            return acb_mat([[src[i, j] for j in indices] for i in range(size)])
        frames.append((columns(kept), columns(omitted)))
    (up, _), (down, e) = frames
    c = adj(step**2)*up
    f = step*e
    projected_f = (eye(19) - e*adj(e))*f
    s = adj(e)*c
    bounds = (arb(".035"), arb(".31"))
    for gram, bound in zip((adj(projected_f)*projected_f, s*adj(s)), bounds, strict=True):
        positive(gram - bound**2*eye(6))
    print("源列及补空间维数", up.ncols(), down.ncols(), e.ncols())
    print("两个最小奇异值的严格下界", *[v.str(5) for v in bounds], flush=True)
    number = 2_306_888_100_090
    _, dst = lib("adjacent-pvm-transport")["spectrum"](lib("three-ray-chiral-boundary"), 6, number)
    high = acb_mat([[dst[i,j] for j in (2,3,14)] for i in range(size)])
    lower_mass_outside(step, source, levels, sectors, high, PROBES)
    assert 2*13 + 3*13 + 4*2 == 73 and 13 + 2*13 + 3*2 == 45


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
