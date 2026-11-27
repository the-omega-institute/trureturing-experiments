"""第五壳共同第四步精确配平完整夸克读数的十五低层见证。"""
from pathlib import Path
from runpy import run_path

from flint import arb


NUMBER = 2_306_888_100_090
QUATERNION = (-11396407304, -13749932891, -73029510396, 997169871114)
DIRECTIONS = ((3, 0, "r"), (14, 9, "r"), (17, 6, "r"), (7, 5, "i"),
              (9, 7, "r"), (10, 9, "i"), (9, 8, "r"), (11, 1, "i"),
              (9, 8, "i"), (15, 7, "r"), (15, 14, "i"), (17, 9, "r"))


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    sites = lib("fifth-flavor-cutoff-chart")["solve"](
        root, "ch39-fifth-flavor-fifteen-low", QUATERNION, DIRECTIONS)
    for i, site in enumerate(sites, 1):
        print("层与截止平方", i, site.str(9), flush=True)
    assert [i+1 for i, s in enumerate(sites) if s > arb("1e-6")] == [1, 3, 4, 15]
    assert all(site < arb("1e-22") for i, site in enumerate(sites, 1) if i not in (1, 3, 4, 15))
    for i, cap in zip((15, 16, 17, 18), ("1.1e-33", "2.4e-34", "5.2e-34", "3.9e-30"), strict=True):
        assert sites[i] < arb(cap)
    weights = next(row[2][2] for row in lib("native-rh-ray-pvm")["shell_weights"](root, NUMBER) if row[0] == 6)
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](root, [NUMBER])
    m = []
    for i, site in enumerate(sites, 1):
        if i in (1, 3, 4, 15):
            m.append(moments(field, heat, site))
        else:
            # 对非负质量有 E_2<=1 与 E_3<=1/2，避免在极低层对区间作除法。
            x = site.upper()
            k_cap = x*sum(abs(a) for _, a, _ in field)/(2*arb.pi())
            v_cap = x*x*sum(abs(g) for g, _, _ in field)/(8*arb.pi())
            m.append((arb(0, k_cap.upper()), arb(0, v_cap.upper())))
    thermal = [sum(c.real*v[j] for c, v in zip(weights, m, strict=True)) for j in range(2)]
    print("牛顿项与真空项", *[v.str(10) for v in thermal], flush=True)
    assert -25 < thermal[0] < -24 and arb("1.6") < thermal[1] < arb("1.7")
    return sites


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
