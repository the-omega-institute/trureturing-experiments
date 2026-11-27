"""两个实际夸克零点上任一单部门 Gram 算子的正谱函数均无法制备六质量。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    model = run_path(str(Path(__file__).with_name("ch39-two-zero-quark-window.py")))
    adj = model["adj"]
    matrices = [model["prepare"](q,arb(x),arb(y)) for q,x,y in (
        (-4,"-0.680061861483031","0.873223965353125"),
        (2,"0.565247837594779","-0.630343671648397"),
    )]
    gram = [adj(a)*a for a in matrices]
    factor = abs(matrices[1].det()/matrices[0].det())**(arb(2)/3)
    boxes = [v+arb(0,w.upper()) for v,w in zip(model["ref"][:6],model["tol"][:6],strict=True)]
    masses = (boxes[:3], boxes[3:])

    def shape(values):
        return sum(x*x for x in values)/(values[0]*values[1]*values[2])**(arb(2)/3)

    required = shape(masses[0])/shape(masses[1])
    assert arb("74.25") < required < arb("86.17")
    for q in (0,1):
        levels,frame = model["eig"](gram[q])
        other = adj(frame)*gram[1-q]*frame
        squared = [x*x for x in masses[q]]
        for assignment in permutations(range(3)):
            ratio = sum(other[i,i].real*squared[j]/levels[i]
                        for i,j in enumerate(assignment))/sum(squared)
            trace_ratio = factor / ratio if q == 0 else factor*ratio
            if q == 0:
                assert trace_ratio < arb("4.884")
            elif assignment[0] == 2:
                assert trace_ratio > 6950
            else:
                assert trace_ratio < arb("5.147")
            print("指定部门、质量排列、归一化迹比上下界",q,assignment,
                  trace_ratio.lower().str(10),trace_ratio.upper().str(10))


if __name__ == "__main__":
    verify()
