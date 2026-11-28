"""用无钟相位的原磁壳固定十二次作用，比较两个中央与具体低编号。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


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
    for q,center,rival,margin in ((-4,2345666,4,305533),(2,964890,11,5620)):
        def score(index):
            return sum(prod((e-r)**2 for r in roots) for e in light(q,index))
        central, other = score(center), score(rival)
        assert central-other > arb(margin)
        print("荷、中央与竞争编号的作用",q,center,central.str(14),rival,other.str(14))


if __name__ == "__main__":
    verify()
