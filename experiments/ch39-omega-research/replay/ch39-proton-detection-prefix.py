"""实际正零点的有限前缀检出固定真空权下的任意非负逐面顶点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def samples(root, count):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))
    shell = load("three-ray-chiral-boundary")
    lib = load("mincut-global-winding")
    target = load("height-independent-flavor-vacuum")
    angles = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    geometry = load("mincut-real-triangle-tomography")
    types, refs, _ = load("fixed-word-proton-faces")["signals"](root)
    ctx.prec = 256
    x = arb("-.01046")
    y = x-angles["reference_angles"](flavor)[0]
    weights = target["minimum_order_weights"](lib,lib["terms"](shell),x,y)
    offset = shell["soft"]["point"]["study"]["offset"]
    rows = []
    for zero in acb.zeta_zeros(1, count):
        h = zero.imag
        delta = arb(23)/25*arb(2)**-20/(h*h+offset)
        phases = shell["geo"]["phases"](h)
        row = []
        for (a,b),label,soft in types:
            loop = geometry["unit"](label)*phases[724]**a*phases[804]**b
            row.append((1-delta*soft)*sum(w*(loop**q).real for q,w in weights.items()))
        rows.append(row)
    return types, refs, rows


def certify(rows, fixed_rows, allowed, witness, support, active, threshold):
    m = len(rows)
    base = [arb(1)] + [arb(0)]*(len(fixed_rows)-1)
    # 前 m-1 枚的精确零由一个方阵的正解给出。
    smaller = fixed_rows + rows[:-1]
    matrix = arb_mat([[row[i] for i in witness] for row in smaller])
    previous = matrix.solve(arb_mat([[v] for v in base+[arb(0)]*(m-1)])).entries()
    assert min(previous) > 0

    # 最优解的活动行编号在 [0,m) 表示 F_i=z，在 [m,2m) 表示 -F_i=z。
    signed = [[value if j < m else -value for value in rows[j % m]] for j in active]
    restricted = [[row[i] for i in support]+[0] for row in fixed_rows]
    restricted += [[row[i] for i in support]+[-1] for row in signed]
    matrix = arb_mat(restricted)
    primal = matrix.solve(arb_mat([[v] for v in base+[arb(0)]*len(active)])).entries()
    p, value = primal[:-1], primal[-1]
    assert min(p) > 0 and abs(value-arb(threshold)) < arb("1e-12")
    for j, row in enumerate(rows):
        fitted = sum(row[i]*mass for i,mass in zip(support,p))
        if j not in active and j+m not in active:
            assert abs(fitted) < value
    dual = matrix.transpose().solve(arb_mat([[0]]*len(support)+[[1]])).entries()
    equality = dual[:len(base)]
    multipliers = [-v for v in dual[len(base):]]
    assert min(multipliers) > 0 and (sum(multipliers)-1).contains(0)
    for i in allowed:
        if i not in support:
            cost = sum(v*row[i] for v,row in zip(multipliers,signed))
            offset = sum(v*row[i] for v,row in zip(equality,fixed_rows))
            assert cost > offset
    assert (equality[0]-value).contains(0)
    print("首个检出前缀、最佳最大读数、前一前缀见证的最小正质量",
          m,value.str(17),min(previous).str(14))


def verify():
    types, refs, rows = samples(Path(__file__).resolve().parent,28)
    fixed_rows = [[arb(1)]*len(types)]
    certify(rows[:17],fixed_rows,range(32),
            (0,6,8,10,11,12,13,14,15,16,17,18,19,20,21,23,26),
            (0,5,7,8,9,11,13,14,17,18,19,20,23,25,26,27),
            (0,1,3,5,6,7,8,9,13,14,15,16,19,21,27,29),".00145502104901")
    # 指定十个实际编号的有理线性组合在每个面轨道上均为正。
    indices = (4,11,13,14,16,19,21,24,25,28)
    multiplier = (11349198,-6194127,-13197936,11293344,7608212,
                  -7707194,14562090,-10909066,-6077534,-11101298)
    scores = [sum(a*rows[i-1][j] for i,a in zip(indices,multiplier)) / sum(map(abs,multiplier))
              for j in range(32)]
    assert min(scores) > arb(".00087193")
    print("十个实际编号的逐点最大模下界",indices,min(scores).str(12))
    allowed = [i for i,key in enumerate(types) if not key[2]]
    assert set(range(32))-set(allowed)=={2,21} and all(refs[i][1]<0 for i in (2,21))
    fixed_rows += [[v[j] for v in refs] for j in (0,2)]
    certify(rows[:16],fixed_rows,allowed,
            (0,5,6,7,8,9,11,12,13,14,15,17,18,19,20,23,26,30),
            (0,6,7,8,9,11,13,14,17,18,19,20,23,25,26,27),
            (0,1,3,6,7,8,9,13,14,15,18,20,26,28),".00063543608986")


if __name__ == "__main__":
    verify()
