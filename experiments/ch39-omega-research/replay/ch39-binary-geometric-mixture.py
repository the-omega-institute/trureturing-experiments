"""一个二值共享变量足以精确实现两夸克中央的十二次校准选择器（至常数）。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_poly, ctx


# 每行是一个共享分支中的下型条件律，整数分子均以一万为分母。
DOWN = (
    ((60,"I",1088), (92,"I",4080), (120,"I",3007), (125,"I",1825)),
    ((1,"I",1758), (37,"I",3330), (68,"I",4908), (76,"D",4)),
)
UP = (
    ((3,"I"), (29,"I"), (31,"I"), (63,"D"), (111,"I"), (125,"D")),
    ((1,"I"), (4,"D"), (60,"I"), (92,"D"), (108,"I"), (124,"D")),
)


def verify():
    replay = Path(__file__).resolve().parent
    shared = run_path(str(replay/"ch39-geometric-selector-mixture.py"))
    selector = run_path(str(replay/"ch39-rational-unique-light-action.py"))
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    ctx.prec = 640
    triple = shared["geometry"](shell, pump)
    factor = shared["squared_factors"]
    target = factor([arb(r) for r in selector["ROOTS"]])
    columns = []
    for up, down in zip(UP, DOWN, strict=True):
        assert all(n > 0 for _,_,n in down) and sum(n for _,_,n in down) == 10000
        mean = sum((arb(n)/10000*factor(triple(2,mask,mode)) for mask,mode,n in down),
                   start=arb_poly([0]))
        columns.extend(mean*factor(triple(-4,mask,mode)) for mask,mode in up)
    weights, minimum, _ = shared["barycentric"](shared["simplex_inverse"](columns), target)
    # 每组六个上型列的权重和为该分支概率；除以该和即得分支内的上型条件律。
    first_mass = sum(weights[i,0] for i in range(len(UP[0])))
    assert 0 < first_mass < 1 and minimum > 0
    # 两组的上下型支撑分别不交，任一侧观察到配方即可恢复分支。
    assert set(UP[0]).isdisjoint(UP[1])
    down_support = [{(m,mode) for m,mode,_ in group} for group in DOWN]
    assert down_support[0].isdisjoint(down_support[1])
    information = -first_mass*first_mass.log()-(1-first_mass)*(1-first_mass).log()
    assert information < arb("0.120756")
    print("两荷配方互信息（自然对数）", information.str(12))
    _, error = shared["rationalize"](weights, columns, target)
    assert error < arb("5.47e-76")
    print("第一分支概率、最小列权、全编号有理化误差上界",
          first_mass.str(12), minimum.str(12), error.str(12))


if __name__ == "__main__":
    verify()
