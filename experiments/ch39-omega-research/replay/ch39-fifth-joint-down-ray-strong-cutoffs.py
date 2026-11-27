"""第五壳六个强截止均低时，原完整味与五式见证的下型射线须移动。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


# 每行给出下型源层、同层组的复整数坐标、荷六重源探针的两个角、
# 要求低截止的荷六层、射线须移动的角下界、固定射线时的最大截止平方下界。
CASES = (
    ((1, 16, 18), ((4425, 14116), (727, 7895), (10733, -17597)),
     "-.05936", "-.41804", (1, 10, 16, 18), ".13240", ".0005838"),
    ((1, 16, 18), ((-19302, 11970), (-3523, 1180), (26747, 5)),
     "3.06958", "-.58206", (16, 18), ".12470", ".0015557"),
    ((1, 16, 18), ((-11724, -1441), (-5023, -2087), (25506, -1821)),
     "-.09920", "-.82683", (18,), ".03029", ".00006870"),
    ((3, 5, 9, 13, 15, 19), ((4993, -4784), (1675, 15751), (9519, 19136),
                             (5842, -2278), (-63, 14195), (17801, -2737)),
     "1.45361", ".12431", (19,), ".02496", ".00008570"),
)


def adj(a):
    return a.conjugate().transpose()


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    center, radius, sectors, source, (m, h, k) = lib("fifth-third-pivot-joint")["verify"](
        root, return_transport=True
    )
    down, _, weights, _, _ = sectors[1]
    _, target = lib("adjacent-pvm-transport")["spectrum"](
        lib("three-ray-chiral-boundary"), 6, 2_306_888_100_090
    )
    for group, numbers, theta, phi, layers, rotation_floor, cutoff_floor in CASES:
        assert tuple(i for i, w in enumerate(weights, 1) if w.overlaps(weights[group[0] - 1])) == group
        basis = acb_mat([[down[i, j - 1] for j in group] for i in range(19)])
        direction = acb_mat([[acb(a, b)] for a, b in numbers])
        e = basis * direction / arb(sum(a*a + b*b for a, b in numbers)).sqrt()
        angle = arb(theta)
        v1, v2 = angle.cos(), angle.sin() * acb(0, arb(phi)).exp()
        probe = source * acb_mat([[v1], [v2]])
        old = center * e
        projected = adj(target) * old
        # 真根的算子误差至多 radius，单位向量的标量读数和投影范数均为 1-Lipschitz。
        incoming = abs((adj(e) * probe)[0, 0]).lower().acos()
        one_step = (abs((adj(e) * old)[0, 0]).lower() - radius).acos()
        projection = sum(abs(projected[i - 1, 0]).lower()**2 for i in layers).sqrt()
        to_layers = (projection - radius).acos()
        distance = (incoming + 2 * one_step + to_layers).upper()
        assert 0 < distance < arb.pi() / 2
        coefficient = len(layers) * (abs(v1)**2 / (h-m) + abs(v2)**2 / (k-m))
        cap = ((arb(".001") - m) * coefficient).sqrt()
        rotation = (cap.acos() - distance) / 3
        fixed_cutoff = (m + distance.cos()**2 / coefficient)**2
        print("荷六低层、下型射线角下界、固定该射线时最大截止平方下界",
              layers, rotation.lower().str(15), fixed_cutoff.lower().str(15), flush=True)
        assert rotation > arb(rotation_floor) and fixed_cutoff > arb(cutoff_floor)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
