"""十二个钟差面的实振幅无法同时抹去两个绕数，以及共同真空的谱下界。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def pair_minimum(a, b):
    # 两个绝对仿射函数的最大值仅在端点、零点和等模交点处需要比较。
    points = [arb(0), arb(1)]
    for left, right in (a, b, (a[0]+b[0], a[1]+b[1]), (a[0]-b[0], a[1]-b[1])):
        assert not (right-left).contains(0)
        at = right/(right-left)
        if 0 < at < 1:
            points.append(at)
        else:
            assert at < 0 or at > 1
    def error(t):
        return abs(t*a[0]+(1-t)*a[1]).max(abs(t*b[0]+(1-t)*b[1]))
    return min(error(t) for t in points)


def verify():
    root = Path(__file__).resolve().parent
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))
    shell = load("three-ray-chiral-boundary")
    formula = load("mincut-global-winding")
    fixed = load("height-independent-flavor-vacuum")
    geometry = load("mincut-real-triangle-tomography")
    angle = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    ctx.prec = 224
    labels = {pair for clock,pair in geometry["enumerate_faces"]().values() if clock == (1,-1)}
    assert labels == {(a,s*b) for a,b in ((13872,23),(-26877,115)) for s in (-1,1)}
    u, v = [geometry["unit"](z) for z in ((13872,23),(-26877,115))]
    rows = [[(u**q).real,(v**q).real] for q in range(1,106)]
    pairs = [(pair_minimum(a,b),i,j) for (i,a),(j,b) in combinations(enumerate(rows,1),2)]
    for selected, answer, floor in ((pairs,(1,18),".00050402"),
                                    ([row for row in pairs if row[2] == 105],(35,105),".00869752")):
        best = min(selected,key=lambda row:float(row[0]))
        assert best[1:] == answer and all(row is best or best[0] < row[0] for row in selected)
        assert best[0] > arb(floor)
        print("最难分开的两个绕数、最优最大实部",answer,best[0].str(16))

    tri = formula["terms"](shell)
    x = -arb(".00815").union(arb(".01280"))
    y = x-angle["reference_angles"](flavor)[0]
    soft_x = [fixed["gradient_columns"](formula,tri,q,x,y)[2] for q in range(1,106)]
    assert soft_x[-1] < -240 and max(abs(v).upper() for v in soft_x[:-1]) < 353
    second = (1-arb(353)/593)/104
    rms = second*arb(".00869752")/arb(2).sqrt()
    assert rms > arb(".0000239")
    print("末圈外至少一项的份额、单位钟差面份额的均方根下界",second.str(14),rms.str(14))


if __name__ == "__main__":
    verify()
