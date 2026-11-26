"""让三个指定实际零点的原始质子因子在同一真空权和同一面概率下同时精确为零。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


NUMBERS = (2_407_782,2_492_462,1_308_440)


def data(root):
    def load(name):
        return run_path(str(root/("ch39-"+name+".py")))
    geo = load("mincut-real-triangle-tomography")
    local = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    fixed = load("height-independent-flavor-vacuum")
    model = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    global_model = load("reference-flavor-global-winding")
    types, _, _ = load("fixed-word-proton-faces")["signals"](root)
    ctx.prec = 224
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    tri = model["terms"](shell)
    weights = fixed["minimum_order_weights"](model,tri,x,y)
    gradient = {q: [arb(1), *fixed["gradient_columns"](model,tri,q,x,y)] for q in range(1,106)}
    responses = []
    for n in NUMBERS:
        h = acb.zeta_zero(n).imag
        clocks = shell["geo"]["phases"](h)
        delta = arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
        row = []
        for (a,b),phase,soft in types:
            z = geo["unit"](phase)*clocks[724]**a*clocks[804]**b
            row.append([arb(0)]+[(1-delta*soft)*(z**q).real for q in range(1,106)])
        responses.append(row)
    return weights, gradient, responses, types, (model,tri,global_model,x,y)


def construct(p, solved, weights, gradient, responses, full=False):
    orders = list(range(1,106)) if full else sorted(set(weights)|set(solved))
    coeff = {q: [*gradient[q], *[sum(v*row[j][q] for j,v in p.items()) for row in responses]] for q in orders}
    fixed = {q: weights.get(q,arb("1e-8")) for q in orders if q not in solved}
    if full:
        fixed = {q: v+arb(0,arb("1e-11")) for q,v in fixed.items()}
    matrix = arb_mat([coeff[q] for q in solved]).transpose()
    assert not matrix.det().contains(0)
    rhs = [arb(int(k==0))-sum(fixed[q]*coeff[q][k] for q in fixed) for k in range(8)]
    # 在中心矩阵预乘后求小修正，避免驻点行与逐面读数行量级相差过大时放宽区间。
    pre = arb_mat([[matrix[i,j].mid() for j in range(8)] for i in range(8)]).inv()
    rhs = arb_mat([[v] for v in rhs])
    center = arb_mat([[z.mid()] for z in (pre*rhs).entries()])
    answer = center+(pre*matrix).solve(pre*(rhs-matrix*center))
    fixed.update(zip(solved,answer.entries()))
    minimum = min(v.lower() for v in fixed.values())
    assert minimum>arb("8.9e-9"), (minimum,[(q,v.str(7)) for q,v in fixed.items() if not v>0])
    for k in range(8):
        actual = sum(fixed[q]*coeff[q][k] for q in orders)
        assert (actual-int(k==0)).contains(0)
    tv=sum(abs(fixed[q]-weights.get(q,0)) for q in orders)/2
    print("参与绕数、最小权、绕数总变差、八阶行列式",len(orders),minimum.str(12),
          tv.str(15),matrix.det().str(12),flush=True)
    return fixed, tv


def vacuum(weights, info):
    model, tri, global_model, x, y = info
    cap = arb(23)/25*arb(2)**-20/14**2
    modes = [(q,[v*c for c in model["coeff"](tri,q,arb(0).union(cap))]) for q,v in weights.items()]
    drift = global_model["global_cover"](model,tri,weights,modes,x,y)
    assert drift < 0
    return drift


def counted_face_floor(gradient, responses, info):
    model, tri, _, x, y = info
    # 全高度共同极小点在单位辐边极限的势能也不超过以下八个有理相位。
    comparisons = [(3157663,293679500,390017), (3255980,-298642388,779346),
                   (2482978,-3967547,21020533), (14006829,-312351443,7293562),
                   (2943189,203091474,1182448), (-14816655,-5324793,7500162),
                   (3191749,-5200950,1722786), (3660583,13718967,2279536)]
    mu = (863157804,-49866606,-86975590)
    direction = (1770684,-42860,769684,228778)
    assert sum(map(abs,mu))==10**9
    columns = []
    for q in range(1,106):
        c = model["coeff"](tri,q,arb(0))
        value = model["energy"](q,c,x,y)
        row = [sum(r[j][q] for j in range(32))/32 for r in responses]
        score = sum(a*b for a,b in zip(mu,row))
        score += sum(a*b for a,b in zip(direction,gradient[q][1:]))
        for u,v,weight in comparisons:
            assert weight>=0
            rival = model["energy"](q,c,arb(u)/10**8,arb(v)/10**8)
            score += weight*(value-rival)
        columns.append(score/10**9)
    bound = min(z.lower() for z in columns)
    assert bound>arb(".007333")
    print("按面原生计数的三编号在所有合格权上的统一读数下界",bound.str(18),flush=True)
    return bound


def verify(root):
    weights,gradient,responses,types,info = data(root)
    info[2]["conjugation_symmetry"](root)
    assert len(types)==32 and types[4]==((0,1),(21097,-103),0)
    sparse14 = {5:arb(20287)/10**5, 27:arb(12443)/10**5, 31:arb(67270)/10**5}
    sparse19 = {5:arb(20078)/10**5, 27:arb(12381)/10**5, 31:arb(67541)/10**5}
    # 全支撑面概率围绕 99.9% 三面分布加 0.1% 原生计数。前 31 个概率独立变，最后一个代回。
    triple = {5:arb(8)/40, 27:arb(5)/40, 31:arb(27)/40}
    spread = {j: (arb(999)/1000*triple.get(j,0)+arb(1)/32000)+arb(0,arb("1e-12")) for j in range(31)}
    spread[31] = 1-sum(spread.values())
    assert min(v.lower() for v in spread.values())>arb(".0000312")
    cases = [
        ("单硬面", {4:arb(1)}, (2,28,30,40,65,67,83,104), False, ".03417"),
        ("原十四圈的近碰撞", sparse14, (1,28,30,65,67,83,104,105), False, ".002711"),
        ("十九圈的近碰撞", sparse19, (1,3,41,45,56,75,79,104), False, ".002559"),
        ("面与绕数双全支撑", spread, (2,3,28,35,38,65,104,105), True, ".02"),
    ]
    for name, probabilities, solved, full, cap in cases:
        assert (sum(probabilities.values())-1).contains(0)
        print("面概率类别、所属轨道数、代回绕数",name,len(probabilities),solved,flush=True)
        obtained, tv = construct(probabilities,solved,weights,gradient,responses,full)
        assert tv<arb(cap)
        drift = vacuum(obtained,info)
        assert drift < -3
    counted_face_floor(gradient,responses,info)


if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
