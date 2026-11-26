"""近真空开窗内作用的精确辨识及有限精度的必要条件。"""
from pathlib import Path
from runpy import run_path
from math import factorial

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    lib = run_path(str(root/"ch39-mincut-global-winding.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 256
    tri = lib["terms"](shell)
    modes = ((1,0), (0,1), (1,-1))
    assert set(tri)=={(0,0), *modes}
    assert sum(len(tri[mode]) for mode in modes)==96
    assert all(n in (0,1) for mode in modes for _,n in tri[mode])
    assert all(n <= 2 for _,n in tri[0,0])
    least = min(max(abs(v).abs_lower() for v in lib["coeff"](tri,q,arb(0)))
                for q in range(1,106))
    assert least > 1
    # 归一，三个径向常数，四个一阶导，以及两组二至十七阶的三个方向。
    constraints = 1+3+4+2*(17-1)*3
    assert constraints==104 < 105
    difference = 2*arb("1e-9")
    angle_radius = arb("9e-6")+arb("2e-7")
    envelope = difference*105*96*(105*angle_radius)**18/factorial(18)
    assert envelope < arb("2e-75")
    print("单位壳最弱绕数的最大谐波振幅下界、线性约束数、权差与作用差上界",
          least.str(13), constraints, difference.str(9), envelope.str(15))


if __name__=="__main__":
    verify()
