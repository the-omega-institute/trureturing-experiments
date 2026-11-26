"""三枚实际零点的同一整数分离证书，在绕数权单纯形与共同真空中的精确总变差半径。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


INDICES = (2_407_782, 2_492_462, 1_308_440)
COEFFICIENTS = (353198, 256731, -390071)
ALTERNATIVE = (34860, 25493, -39647)


def data(root):
    def load(name):
        return run_path(str(root / ("ch39-"+name+".py")))
    types, _, _ = load("fixed-word-proton-faces")["signals"](root)
    geometry = load("mincut-real-triangle-tomography")
    model = load("mincut-global-winding")
    local = load("mincut-reference-ckm-action")
    shell = load("three-ray-chiral-boundary")
    flavor = load("chiral-neutral-flavor-phase")
    fixed = load("height-independent-flavor-vacuum")
    global_model = load("reference-flavor-global-winding")
    ctx.prec = 224
    x = arb("-.01046")
    y = x-local["reference_angles"](flavor)[0]
    tri = model["terms"](shell)
    weights = fixed["minimum_order_weights"](model, tri, x, y)
    gradients = {q: [arb(1), *fixed["gradient_columns"](model,tri,q,x,y)] for q in range(1,106)}
    scores, alternate = ([[arb(0)]*106 for _ in types] for _ in range(2))
    for n, integer, second in zip(INDICES, COEFFICIENTS, ALTERNATIVE):
        h = acb.zeta_zero(n).imag
        delta = arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
        phases = shell["geo"]["phases"](h)
        for j, ((a,b), label, soft) in enumerate(types):
            loop = geometry["unit"](label)*phases[724]**a*phases[804]**b
            factor = arb(integer)*(1-delta*soft)/sum(map(abs,COEFFICIENTS))
            other = arb(second)*(1-delta*soft)/sum(map(abs,ALTERNATIVE))
            for q in range(1,106):
                amplitude = (loop**q).real
                scores[j][q] += factor*amplitude
                alternate[j][q] += other*amplitude
    return weights, gradients, scores, alternate, (model, tri, global_model, x, y)


def transfers(row, weights, gradients, add, remove, stationary):
    support = list(add)+list(remove)
    signs = [1]*len(add)+[-1]*len(remove)
    def column(q):
        return (gradients[q] if stationary else [arb(1)]) + [row[q]]
    mat = arb_mat([[sign*v for v in column(q)] for q,sign in zip(support,signs)]).transpose()
    center = sum(weights[q]*row[q] for q in weights)
    rhs = [arb(0)]*(mat.nrows()-1)+[-center]
    amounts = mat.solve(arb_mat([[v] for v in rhs])).entries()
    assert min(amounts)>0
    final = weights.copy()
    for q, sign, amount in zip(support,signs,amounts):
        final[q]=final.get(q,arb(0))+sign*amount
    assert min(final.values())>0
    radius = sum(amounts[len(add):])
    dual = mat.transpose().solve(arb_mat([[arb(sign<0)] for sign in signs])).entries()
    return final, radius, amounts, dual, column


def vacuum(weights, info):
    model, triangles, global_model, x, y = info
    radial_cap = arb(23)/25*arb(2)**-20/14**2
    deficit = arb(0).union(radial_cap)
    modes = [(q, [v*c for c in model["coeff"](triangles, q, deficit)]) for q,v in weights.items()]
    drift = global_model["global_cover"](model,triangles,weights,modes,x,y)
    assert drift < 0
    return drift


def verify_boundary(weights, gradients, scores, alternate, info, allowed, stationary):
    # 正质量从原支持中移走，可以加到允许的任何绕数上。各列的最坏下降斜率给闭 TV 球上的余量。
    extended = len(allowed)==105
    special = {31: ((60 if extended else 67,), (1,))}
    if stationary:
        if extended:
            special = {31: ((41,43,75,77,104),(1,)),
                       5: ((10,28,45),(40,83,105)),
                       6: ((37,98,101,104),(4,35))}
        else:
            special = {31: ((30,67,83,105),(1,104))}
    bounds = []
    costs = []
    centers = []
    for j,row in enumerate(scores):
        center = sum(w*row[q] for q,w in weights.items())
        centers.append(center)
        assert center > arb(".00257774")
        if j in special:
            add, remove = special[j]
            boundary, radius, amounts, dual, column = transfers(row,weights,gradients,add,remove,stationary)
            # 质量增加列的成本为零，移走列的成本为一。基列以方程恒等取等，其余严格验对偶可行。
            def reduced(q):
                return sum(a*b for a,b in zip(dual,column(q)))
            assert all(reduced(q)<0 for q in allowed if q not in add)
            assert all(reduced(q)>-1 for q in weights if q not in remove)
            assert dual[-1]<0 and (radius+center*dual[-1]).contains(0)
            cost = -1/dual[-1]
            if j==31:
                hit = boundary
                winning = radius
                print("最先触零的轨道、质量增加与移走绕数、转移量",
                      j,add,remove,[v.str(12) for v in amounts],flush=True)
        else:
            # 忽略驻点约束仍可排除其余列，这给其可下降速率的上界。
            cost = max(row[q].upper() for q in weights)-min(row[q].lower() for q in allowed)
            assert cost>0
            radius = center/cost
        bounds.append(radius)
        costs.append(cost)
    assert min(bounds[j] for j in range(32) if j!=31)>winning
    assert (sum(hit[q]*scores[31][q] for q in hit)).contains(0)
    assert (sum(hit.values())-1).contains(0)
    if stationary:
        for k in range(1,5):
            assert sum(hit[q]*gradients[q][k] for q in hit).contains(0)
        drift = vacuum(hit,info)
        assert drift<arb("-3.50")
        alternative_floor = min(sum(hit[q]*row[q] for q in hit).lower() for row in alternate)
        print("同一三编号在边界权上改用另一整数乘子的逐面余量",alternative_floor.str(16),flush=True)
        assert alternative_floor>arb(".0004")
    r=arb(".0015")
    at_r=min((center-r*cost).lower() for center,cost in zip(centers,costs))
    print("允许绕数、共同真空约束、锐 TV 半径、最坏轨道中心余量、TV ≤ 0.0015 时统一检出余量",
          len(allowed),stationary,winning.str(22),centers[31].str(16),at_r.str(16),flush=True)
    assert at_r>0
    return winning, at_r


def full_support(weights,gradients,scores):
    # 复用百维共同全局真空族定义的精确参数盒，估其到十四项中心的统一 TV 距离。
    solved = (30,35,38,40,105)
    free = [q for q in range(1,106) if q not in solved]
    coeff = arb_mat([gradients[q] for q in solved]).transpose()
    box = {q: weights.get(q,arb("1e-8"))+arb(0,arb("1e-9")) for q in free}
    rhs = [arb(int(i==0))-sum(box[q]*gradients[q][i] for q in free) for i in range(5)]
    box.update(zip(solved,coeff.solve(arb_mat([[v] for v in rhs])).entries()))
    assert min(box.values())>arb("8.9e-9")
    tv=sum(abs(v-weights.get(q,0)) for q,v in box.items()).upper()/2
    margin=min(sum(row[q]*v for q,v in box.items()).lower() for row in scores)
    print("已证百维全正共同真空盒到原十四权的 TV 上界及三枚证书的统一余量",tv.str(18),margin.str(18),flush=True)
    assert tv<arb(".0001") and margin>arb(".0025")


def verify(root):
    weights, gradients, scores, alternate, info = data(root)
    info[2]["conjugation_symmetry"](root)
    for stationary in (False,True):
        for allowed in (sorted(weights),range(1,106)):
            radius, margin = verify_boundary(weights,gradients,scores,alternate,info,allowed,stationary)
            expected = { (False,14): ".00193342831275027", (False,105): ".00155814271316563",
                         (True,14): ".00212037598012562", (True,105): ".00191778780666884" }
            assert abs(radius-arb(expected[stationary,len(allowed)]))<arb("1e-16")
    full_support(weights,gradients,scores)


if __name__=="__main__":
    root = Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
