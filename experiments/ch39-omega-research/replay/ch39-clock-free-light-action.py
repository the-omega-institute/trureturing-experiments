"""用无钟相位的原磁壳固定十二次作用，比较两个中央与具体低编号。"""
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def verify():
    replay = Path(__file__).resolve().parent
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    lowest = run_path(str(replay/"ch39-global-light-action-degree.py"))["lowest"]
    read = run_path(str(replay/"ch39-light-spectrum-zero-selector.py"))
    light = read["light_spectra"](replay)
    ctx.prec = 256
    blocks = list(pump["rotation_bases"](shell))
    adj = pump["adjoint"]
    roots = []
    for charge in (-4, 2):
        base, edges = pump["clock_matrices"](shell,charge)
        link = edges[0]+edges[1]
        flat = base+link+adj(link)
        roots.extend(lowest(adj(b)*flat*b) for b in blocks)
    # 任取单位辐边参考角，三个字符能量不能全落入目标的六个根盒。
    inverse = run_path(str(replay/"ch39-rational-unique-light-action.py"))
    target_radii = arb("1e-23")
    target_boxes = [arb(root)+arb(0,target_radii.upper()) for root in inverse["ROOTS"]]
    wide = arb(-1).union(arb(1))
    zero_radial = acb_mat(19,19)
    for charge in (-4,2):
        atlas = inverse["determinant_rows"](pump,shell,zero_radial,blocks,charge,target_boxes)
        for indices in product(range(6),repeat=3):
            rows = [block[i] for block,i in zip(atlas,indices,strict=True)]
            p,_,_ = inverse["eliminate"](rows,[wide]*3)
            assert not p(0).contains(0)
    print("任意两个单位辐边参考点，各有一条能量避开全部六个有理根盒，半宽",target_radii)
    for q,center,rival,margin in ((-4,2345666,4,305533),(2,964890,11,5620)):
        def score(index):
            return sum(prod((e-r)**2 for r in roots) for e in light(q,index))
        central, other = score(center), score(rival)
        assert central-other > arb(margin)
        print("荷、中央与竞争编号的作用",q,center,central.str(14),rival,other.str(14))


if __name__ == "__main__":
    verify()
