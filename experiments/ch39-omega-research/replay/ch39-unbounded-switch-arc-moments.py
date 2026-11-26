"""四个非中心零点给出可共同平分任意有限批单弧一阶矩的权重。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

fixed = run_path(str(Path(__file__).with_name("ch39-fixed-schedule-identification.py")))
labels = (1, 28814, 97370, 1475240)


def verify():
    ctx.prec = 160
    projection = fixed["flavor"]["arcs"][32]
    complement = fixed["flavor"]["identity"]-projection
    # rank-one 交叉块的 (0,0) 坐标非退化，因而一个复数已确定整块。
    assert 0 < projection[0,0].real < 1
    center_height = acb.zeta_zero(11431).imag
    offset = arb(1)/4+30+arb(2)/15
    reference = fixed["native"]["plane_response"](center_height)
    rows = []
    for label in labels:
        height = acb.zeta_zero(label).imag
        gauge = (center_height**2+offset)/(height**2+offset)-1
        difference = fixed["native"]["plane_response"](height)-reference
        component = (projection*difference*complement)[0,0]
        rows.append([1,gauge,component.real,component.imag])
    features = arb_mat(rows).transpose()
    assert 0 not in features.det()
    weights = features.solve(arb_mat([[1],[0],[0],[0]]))
    assert all(v>0 for v in weights.entries())
    print("Labels",labels)
    print("Positive exact-solution weights",*(weights[j,0].str(17) for j in range(4)))


if __name__ == "__main__":
    verify()
