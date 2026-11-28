"""保留每条线的总实权，在四个线对上分别允许两个正份额截止。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent

# 每行对应第一部门的一条线，每列对应第二部门的一条线。
# 以十二进制数字 k 表示插值 k/11。四条分拆线在下方覆盖。
LABELS = (
    "b000010000010001bba", "2bbbb09b8870744b434", "2bbaa0ba878b652b023",
    "1bab9094600b030b000", "1ba89075310b003b000", "1bb77094530b203b210",
    "0000070bbbb0bbb0bbb", "17847b20000b000b000", "19643b00000b000a001",
    "000ba0babbb0bbb0bbb", "13411b00010a0008000", "17604b10000b000b000",
    "18405b30001b0019010", "17713b10200b000a002", "12500b20000a1009020",
    "13600b000008100b000", "17501b10000900050b2", "0bbbb0bbba90bba0302",
    "15210b11010b000657b",
)
SPLITS = {0: (0, 11), 16: (0, 11), 17: (10, 11), 338: (0, 11)}


def verify():
    previous = run_path(str(ROOT / "ch39-two-copy-gram-cutoff-two-exceptions.py"))
    edges, rows, heat, _ = previous["model_inputs"]()
    lib = run_path(str(ROOT / "ch39-chiral-minimal-neutral-gravity.py"))
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    gauge = run_path(str(ROOT / "ch39-signed-gravity-unification.py"))
    moments = run_path(str(ROOT / "ch39-signed-cutoff-gravity.py"))["moments"]
    fields = gauge["gauge_weights"](model)
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    logscale = (arb(10**16)/model["mz"]).log()
    baseline = arb_mat([[0], [0], *[[v-b*logscale] for v,b in zip(initial, beta, strict=True)]])
    assert len(LABELS) == 19 and all(len(row) == 19 for row in LABELS)
    directions = []
    for i, ((a, b, weight), digit) in enumerate(zip(edges, "".join(LABELS), strict=True)):
        def response(label):
            x = arb(7)/1000*((11-label)*a+label*b)/11
            assert 0 < x <= 1
            data = gauge["threshold_column"](rows, moments, heat, fields, x)[1:]
            return arb_mat([[weight*value] for value in data])

        if i in SPLITS:
            left, right = SPLITS[i]
            base = response(right)
            directions.append(response(left)-base)
        else:
            base = response(int(digit, 12))
        baseline += base
    # 行顺序 K,W,三支逆耦合。最后的两个方程规定三支在同一能标相等。
    projection = arb_mat([[1,0,0,0,0], [0,1,0,0,0], [0,0,1,-1,0], [0,0,1,0,-1]])
    rhs = arb_mat([[1],[arb("1e-122")],[0],[0]]) - projection*baseline
    delta = arb_mat([v.entries() for v in directions]).transpose()
    coefficients = projection*delta
    assert not coefficients.det().contains(0)
    fractions = coefficients.solve(rhs)
    assert all(0 < f < 1 for f in fractions.entries())
    readings = baseline+delta*fractions
    common = readings[4,0]
    assert beta[2] < beta[1] < 0 < beta[0] and logscale > 0
    minimum = common+beta[2]*logscale
    assert arb("37.56") < common < arb("37.57") and minimum > arb("1.54")
    print("四条分拆线在较低位置的份额", *[v.str(14) for v in fractions.entries()])
    print("共同逆耦合、M_Z 到交汇能标的全程最小值", common.str(16), minimum.str(16))


if __name__ == "__main__":
    verify()
