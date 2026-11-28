"""同一有理十四次可加作用在两个味部门及中性壳上分别唯一选择校准编号。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    replay = Path(__file__).resolve().parent
    model = run_path(str(replay/"ch39-rational-unique-light-action.py"))
    neutral = run_path(str(replay/"ch39-quark-action-neutral-height.py"))
    neutral["verify"]()  # 非平凡中性字符能量在所有实际高度严格递增。
    read = run_path(str(replay/"ch39-light-spectrum-zero-selector.py"))
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    lowest = run_path(str(replay/"ch39-global-light-action-degree.py"))["lowest"]
    light = read["light_spectra"](replay)
    ctx.prec = 256
    adj = pump["adjoint"]
    block = list(pump["rotation_bases"](shell))[1]

    def energy(n):
        matrix = shell["laplacian"](0,acb.zeta_zero(n).imag)
        return lowest(adj(block)*matrix*block)

    lower, middle, upper = [energy(n) for n in (11430,11431,11432)]
    neutral_root = arb("8.2468623769185075101687162065427022552841417343")
    assert abs(middle-neutral_root) < arb("1e-45")
    assert middle-lower > arb("6.56e-19") and upper-middle > arb("3.53e-19")
    roots = [arb(s) for s in model["ROOTS"]]+[neutral_root]
    # 先前全高度区间以这个较宽的区间包住中性非平凡最低能量。
    neutral_band = arb("8.24686").union(arb("8.24687"))
    assert all(abs(neutral_band-r) > arb(".20") for r in roots[:-1])
    assert min(neutral_root-lower,upper-neutral_root) > arb("3.52e-19")
    h = lambda x: prod((x-r)**2 for r in roots)
    neutral_gap = 2*((arb("3.52e-19")**2 * arb(".20")**12)-h(middle))
    assert neutral_gap > arb("1e-45")
    print("中性部门扣除共同常数后的严格作用差下界",neutral_gap.str(12))

    eps = arb("1e-32")
    model["verify_number_boxes"](pump,shell,model["CENTERS"],[r+arb(0,eps.upper()) for r in roots])
    low = eps**2 * min(
        prod((abs(x-y)/2)**2 for j,y in enumerate(roots) if j!=i).lower()
        for i,x in enumerate(roots)
    )
    for q,n in model["CENTERS"]:
        gap = low - sum(h(e) for e in light(q,n))
        assert gap > arb("1e-69")
        print("味部门、所有其他编号的严格作用差下界",q,gap.str(12))


if __name__ == "__main__":
    verify()
