"""两枚不同的实际零点在必要高度窗中同时给出零费用，且能被其谱作用唯一选中。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb


PAIR = (2345467,964715)
CLASSES = (
    (1,1,5,1,0,"I","D"), (2,0,1,0,0,"D","I"), (3,0,1,0,0,"D","I"),
    (4,0,1,1,0,"D","I"), (4,0,4,1,4,"I","D"), (4,1,5,1,4,"I","D"),
    (5,0,5,1,5,"D","D"), (5,1,6,0,5,"D","D"), (6,0,4,1,4,"D","D"),
    (6,0,5,0,5,"I","D"), (6,0,5,1,5,"I","D"), (6,0,6,1,6,"I","D"),
)


def verify():
    replay = Path(__file__).resolve().parent
    setup = run_path(str(replay/"ch39-minimum-orbit-specific-probability.py"))["setup"]
    contract = run_path(str(replay/"ch39-common-zero-orbit-cost.py"))["unique_gap"]
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    optimal = setup()
    offsets = (PAIR[0]-2345666,PAIR[1]-964890)
    assert optimal(offsets,[("class",v) for v in CLASSES]).is_zero()
    u,d = (acb.zeta_zero(n).imag for n in PAIR)
    assert u-d > arb("728522.42") and arb("2.25390") < u/d < arb("2.26860")
    gap = contract(replay,shell,PAIR)
    assert gap > arb("6.13e-69")
    print("两个实际编号、高度差、高度比与其他编号作用差下界",
          PAIR,(u-d).str(12),(u/d).str(12),gap.str(12))


if __name__ == "__main__":
    verify()
