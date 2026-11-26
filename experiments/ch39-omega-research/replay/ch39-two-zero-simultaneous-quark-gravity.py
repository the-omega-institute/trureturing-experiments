"""在同一复三代空间和同一正有理读数下，同时核验两个相邻零点的味与引力。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify(root):
    frame,fermions = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))["verify"](root)
    cutoffs = run_path(str(root / "ch39-adjacent-ordered-charge-six-cutoffs.py"))["verify"](root)
    ordered = run_path(str(root / "ch39-bounded-common-rational-amplitude.py"))["ordered"]
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 320
    n = 2_306_888_100_085
    data = []
    for number,(levels,positions) in zip((n,n+1),cutoffs[1:],strict=True):
        height = acb.zeta_zero(number).imag
        assert all(x.sqrt()*arb("1.22089e19")>arb("91.1876") for x in positions)
        data.extend(zip(levels,[x.sqrt() for x in positions],strict=True))
        for charge in (-4,2):
            data.extend(zip(ordered(shell,charge,height),fermions[number,charge],strict=True))
    data.sort(key=lambda row:row[0].mid())
    assert len(data)==114 and all(x<y for (x,w),(y,v) in zip(data,data[1:]))
    floor = arb("91.1876") / arb("1.22089e19")
    assert all(floor<w<arb(".948") for _,w in data)
    min_gap = min((y-x).lower() for (x,w),(y,v) in zip(data,data[1:]))
    print("两编号、三荷 114 个真实节点的最小间距",min_gap.str(13),flush=True)
    left_cut = cutoffs[0][1][0].sqrt()*arb("1.22089e19")
    print("左邻三段解的最低截止（GeV）",left_cut.str(13),flush=True)
    return frame, data


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
