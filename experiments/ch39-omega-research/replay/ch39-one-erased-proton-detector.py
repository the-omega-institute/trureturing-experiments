"""在固定十四项真空权下，认证丢失任一读数后的四枚实际零点检出名单。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


# 每行按剩余三个编号的原顺序给整数对偶，以及最坏分布的三个面轨道。
INDICES = (7_454_773, 1_290_431_208, 416_569_899, 84_227_309)
ERASURES = (
    ((394208, -569384, -36408), (18, 25, 16)),
    ((395939, -39038, -565024), (19, 24, 17)),
    ((404935, -24810, -570256), (19, 24, 7)),
    ((-23434, 403546, -573019), (25, 26, 16)),
)


def read_faces(root):
    def load(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    types, _, _ = load("fixed-word-proton-faces")["signals"](root)
    geometry = load("mincut-real-triangle-tomography")
    model = load("mincut-global-winding")
    angles = load("mincut-reference-ckm-action")
    shell = load("three-ray-chiral-boundary")
    flavor = load("chiral-neutral-flavor-phase")
    fixed = load("height-independent-flavor-vacuum")
    ctx.prec = 256
    x = arb("-.01046")
    y = x - angles["reference_angles"](flavor)[0]
    weights = fixed["minimum_order_weights"](model, model["terms"](shell), x, y)
    assert len(weights) == 14 and min(weights.values()) > 0 and (sum(weights.values())-1).contains(0)
    rows = []
    for n in INDICES:
        height = acb.zeta_zero(n).imag
        delta = arb(23)/25 * arb(2)**-20 / (height**2 + shell["soft"]["point"]["study"]["offset"])
        phases = shell["geo"]["phases"](height)
        row = []
        for (a, b), pair, soft in types:
            loop = geometry["unit"](pair) * phases[724]**a * phases[804]**b
            row.append((1-delta*soft) * sum(w*(loop**q).real for q, w in weights.items()))
        assert len(row) == 32
        rows.append(row)
        print("实际零点、其高度及双钟主相角", n, height.str(20),
              phases[724].arg().str(17), phases[804].arg().str(17))
    assert acb.zeta_zero(1).imag > 14 and shell["soft"]["point"]["study"]["offset"] > 0
    return rows, sum(q*w for q, w in weights.items())


def certify(triple, integers, support):
    norm = sum(map(abs, integers))
    floor = min((sum(a*row[j] for a, row in zip(integers, triple))/norm).lower()
                for j in range(32))
    signed = [[(1 if a > 0 else -1)*v for v in row] for a, row in zip(integers, triple)]
    matrix = arb_mat([[row[k] for k in support]+[-1] for row in signed]+[[1,1,1,0]])
    primal = matrix.solve(arb_mat([[0],[0],[0],[1]])).entries()
    dual = matrix.transpose().solve(arb_mat([[0],[0],[0],[-1]])).entries()
    assert min(primal) > 0 and min(dual[:3]) > 0 and (sum(dual[:3])-1).contains(0)
    assert (primal[-1]+dual[-1]).contains(0)
    for j in range(32):
        if j not in support:
            assert sum(a*row[j] for a, row in zip(dual, signed)) > primal[-1]
    return floor, primal[-1], primal[:3]


def verify(root):
    rows, slope = read_faces(root)
    floors, optima = [], []
    for j, (integers, support) in enumerate(ERASURES):
        floor, optimum, probabilities = certify(rows[:j]+rows[j+1:], integers, support)
        floors.append(floor)
        optima.append(optimum)
        print("遗漏编号、整数证书下界、锐检出值、最坏轨道及组概率", INDICES[j],
              floor.str(17), optimum.str(19), support, [p.str(16) for p in probabilities])
    assert all(optima[2] < v for j, v in enumerate(optima) if j != 2)
    assert arb(".0019807") < optima[2] < arb(".0019809")
    cap = arb(23)/25 * arb(2)**-20/14**2
    robust = min(floors) - 2*arb("1e-5")*slope - cap
    assert robust > arb(".00095900")
    print("任失一枚的锐检出值、各角容许 1e-5 弧度后的统一下界",
          optima[2].str(20), robust.str(16))


if __name__ == "__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
