"""在右外壳上量化六质量和完整 CKM 中心所需的上型四站谱形变。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def reflection_spectrum(root):
    model = run_path(str(root / "ch39-rank-two-generator-up-quark.py"))
    training = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    transport = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    shells = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    orbits = run_path(str(root / "ch39-up-quark-arbitrary-site-phase-floor.py"))["ORBITS"]
    n = model["N"]
    _, responses, _ = training["verify_nodes"](
        root, ((n + 1, ".0294"),), ".1", (), ((6, 7, 10, 13), (6, 7, 10, 12))
    )
    ctx.prec = 384
    labels = responses[n + 1, -4]
    diameter_bound = max(x.upper() for x in labels) - min(x.lower() for x in labels)
    reflection = model["eye"](19) - 2 * model["projector"](orbits)
    _, source = transport["spectrum"](shells, -4, n + 1)
    _, target = transport["spectrum"](shells, -4, n + 2)
    weights = model["read"](source, target, labels, reflection)
    sites = (shells["soft"]["pole"], *shells["boundary_sites"]())
    basis = acb_mat([[target[i,j] for j in range(19)] for i in sites])
    adj = model["adj"]
    raw = acb_mat([[basis[i,j] * weights[j] for j in range(19)] for i in range(4)]) * adj(basis)
    sqrt3 = arb(3).sqrt()
    omega = acb(-arb(1)/2, sqrt3/2)
    characters = []
    for phase in (omega, omega.conjugate()):
        vector = acb_mat([[0], [1], [phase], [phase**2]]) * (1/sqrt3)
        value = (adj(vector) * raw * vector)[0,0].real
        assert (raw * vector - vector * value).contains(acb_mat(4,1))
        characters.append(value)
    trivial = acb_mat([[1,0], [0,1/sqrt3], [0,1/sqrt3], [0,1/sqrt3]])
    block = adj(trivial) * raw * trivial
    center = (block[0,0].real + block[1,1].real)/2
    half_gap = (((block[0,0].real - block[1,1].real)/2)**2 + abs(block[0,1])**2).sqrt()
    return (*characters, center - half_gap, center + half_gap), diameter_bound


def probability_box(reference, radius, masses):
    error = arb(0, radius.upper())
    e, f, low, high = [z + error for z in reference]
    u, c, t = masses
    assert u < e < c < high < t < f
    # 由正性和压缩交错，介于零与最轻上型质量之间的只能是 low。
    low = low.intersection(arb(0).union(u))
    weights = []
    for i, m in enumerate(masses):
        derivative = prod(m - y for j, y in enumerate(masses) if j != i)
        # 分开除以已与零分离的因子，以减少宽区间乘积带来的精度损耗。
        ratio = abs(((m - low)*(m - high)) / (m - e) / (m - f) / derivative)
        assert ratio.is_finite()
        weights.append(ratio)
    lower = [max(arb(0), z.lower()) for z in weights]
    upper = [z.upper() for z in weights]
    probabilities = []
    for i in range(3):
        others = [j for j in range(3) if i != j]
        left = lower[i] / (lower[i] + sum(upper[j] for j in others))
        right = upper[i] / (upper[i] + sum(lower[j] for j in others))
        assert left.is_finite() and right.is_finite()
        probabilities.append(left.union(right))
    return probabilities


def exclude_chart(lib, geom, mixing, base, other, probabilities):
    stack = [(arb(-1), arb(1), (0,1,2), 0)]
    count = depth = 0
    least = arb(1)
    while stack:
        lo, hi, ids, level = stack.pop()
        count += 1
        depth = max(depth, level)
        lines = lib["lines"](base, other, lo, hi, ids, geom)
        keep = []
        for i in ids:
            if i in lines:
                observed = [abs(z)**2 for z in (mixing * lines[i]).entries()]
                gap = max(abs(p - q).lower() for p, q in zip(observed, probabilities, strict=True))
                if gap > 0:
                    least = min(least, gap)
                    continue
            keep.append(i)
        if keep:
            assert level < 24, (keep, lo, hi)
            mid = (lo + hi)/2
            stack.extend(((lo, mid, keep, level + 1), (mid, hi, keep, level + 1)))
    return count, depth, least


def verify(root):
    reference, diameter_bound = reflection_spectrum(root)
    pencil = run_path(str(root / "ch39-triple-geometry-pencil.py"))
    lib, geom, mixing, masses, a, b, _ = pencil["setup"](root)
    radius = arb(".00118")
    probabilities = probability_box(reference, radius, masses[0])
    print("反射处的两字符与平凡块升序谱", [z.str(15) for z in reference])
    for z in probabilities:
        print("上型质量本征线概率区间", z.lower().str(12), z.upper().str(12))
    for name, base, other in (("A+tB", a, b), ("B+tA", b, a)):
        count, depth, least = exclude_chart(lib, geom, mixing, base, other, probabilities)
        print(name, "累计区间、最大深度、概率距离下界", count, depth, least.str(12))
    transport_bound = radius / diameter_bound
    assert transport_bound > arb(".00124531")
    print("谱盒半径、源谱直径上界、去全局相位的输运距离下界",
          radius.str(12), diameter_bound.str(15), transport_bound.str(15))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
