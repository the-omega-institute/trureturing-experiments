"""实际三源经保持三循环的整系数滤波后，检验夸克正扩张和所需的最小形变。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def odds(length):
    return length / (1-length**2).sqrt()


def verify(root):
    free = run_path(str(root / "ch39-free-four-site-positive-extensions.py"))
    mixing = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    gram = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    triplets = run_path(str(root / "ch39-four-gram-triplet-ckm.py"))
    projector = run_path(str(root / "ch39-three-gram-source-flavor.py"))["projections"]
    ctx.prec = 320
    angles = [arb(value) for value in (".22704", ".04275", ".00377", "1.139")]
    target = mixing["standard_matrix"](*[v.sin() for v in angles[:3]], angles[3])
    table = gram["gram_cholesky"](root)
    upper = arb(".0332714")
    lengths = []
    for omit in range(4):
        n = [v.real for v in triplets["normal"](table, omit)]
        fixed = (n[0]**2 + sum(n[1:])**2/3).sqrt()
        lengths.append(fixed)
        required = odds(fixed) / odds(upper)
        print("省略源、可逆等变输运的条件数必要下界", omit, required.lower().str(10))
        assert required > arb("18.27")

    # S=I+10(I+C+C²) 在固定平面上乘 31，在转动平面上乘 1。
    # 作用于三源后，法向随逆转置变换。直接从变换后的三行重算作为核对。
    rows = [[31*row[0], *[x+10*sum(row[1:]) for x in row[1:]]] for row in table]
    n = [v.real for v in triplets["normal"](rows, 2)]
    length = (n[0]**2 + sum(n[1:])**2/3).sqrt()
    old = lengths[2]
    assert length.overlaps(old/(old**2+31**2*(1-old**2)).sqrt())
    assert arb(".02163137") < length < arb(".02163138")
    print("整数滤波后实法向的固定分量", length.str(14))
    free["witness"](target, 2, length, projector, down_weight=arb("1.016334"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
