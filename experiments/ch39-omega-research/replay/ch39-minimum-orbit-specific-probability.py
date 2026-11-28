"""为精确复现中央校准多项式，求必须容许区分轨道身份的最小概率。"""
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly, ctx


def setup():
    replay = Path(__file__).resolve().parent
    def load(name):
        return run_path(str(replay/f"ch39-{name}.py"))

    mixture = load("geometric-selector-mixture")
    grouped = load("geometric-exchangeable-priors")
    reference = load("rational-unique-light-action")
    shell, pump = load("three-ray-chiral-boundary"), load("native-light-projector-pump")
    light = load("light-spectrum-zero-selector")["light_spectra"](replay)
    ctx.prec = 384
    triple = mixture["geometry"](shell,pump)
    factor = mixture["squared_factors"]
    labels = list(product(grouped["MASKS"],grouped["MODES"]))
    up, down = [{(m,b):factor(triple(q,m,b)) for m,b in labels} for q in (-4,2)]
    classes = grouped["geometry_classes"]()
    index = {pair:i for i,pair in enumerate(labels)}
    coeffs = [arb_mat([[poly[label][k] for k in range(7)] for label in labels])
              for poly in (up,down)]

    def optimal(offsets, support):
        recipes, cost = [], []
        for tag, values in support:
            if tag == "class":
                u,pu,d,pd,overlap,b,c = values
                pairs = classes[u,pu,d,pd,overlap]
                recipes.append(sum((up[i,b]*down[j,c] for i,j in pairs),arb_poly([0]))
                               * (arb(1)/len(pairs)))
                cost.append(0)
            else:
                assert tag == "pair"
                i,b,j,c = values
                recipes.append(up[i,b]*down[j,c])
                cost.append(1)
        du, dd = offsets
        levels = light(-4,2345666+du)+light(2,964890+dd)
        ns = [(x*10**46+arb(1)/2).floor().unique_fmpz() for x in levels]
        assert all(n is not None for n in ns)
        roots = [arb(n)/10**46 for n in ns]
        if offsets == (0,0):
            assert [int(n) for n in ns] == [Fraction(t)*10**46 for t in reference["ROOTS"]]
        target = factor(roots)
        inverse = mixture["simplex_inverse"](recipes)
        weights, _, _ = mixture["barycentric"](inverse,target)
        upper = sum(c*weights[i,0] for i,c in enumerate(cost))

        # 对偶仿射函数在这十二列上精确取各自的费用。计算所有顶点上的
        # 最大超额并下移常数项，即便其他列重复或紧贴该面也得到可行下界。
        dual = arb_mat([cost])*inverse
        kernel = arb_mat([[dual[0,i+j-1] if 1 <= i+j <= 11 else 0
                           for j in range(7)] for i in range(7)])
        a,b = coeffs
        scores = a*kernel*b.transpose()
        bias = dual[0,11]
        excess = max(arb(0), (max(v.upper() for v in scores.entries())+bias-1).upper())
        for pairs in classes.values():
            for bu,bd in grouped["MODE_PAIRS"]:
                avg = sum(scores[index[i,bu],index[j,bd]] for i,j in pairs)/len(pairs)
                excess = max(excess,(avg+bias).upper())
        # 同一原始与对偶基给出相等目标值；减去超额后对全凸包成立。
        lower = upper-excess
        assert min(weights[i,0].lower() for i in range(12)) > 0
        assert excess < arb("1e-70")
        return lower.lower().union(upper.upper())
    return optimal


def verify():
    optimal = setup()
    data = json.loads(Path(__file__).with_name("ch39-nine-orbit-cost-bases.json").read_text())
    costs = {tuple(item["offsets"]):optimal(tuple(item["offsets"]), item["support"])
             for item in data}
    assert set(costs) == set(product((-1,0,1),repeat=2))
    center = costs[0,0]
    assert arb("0.00016543237") < center < arb("0.00016543240")
    assert all(cost > 2*center for k,cost in costs.items() if k != (0,0))
    saved = Path(__file__).with_name("ch39-holdout-orbit-cost-bases.json")
    held = {tuple(item["offsets"]):optimal(tuple(item["offsets"]), item["support"])
            for item in json.loads(saved.read_text())}
    expected = {pair for pair in product(range(-2,3),repeat=2) if max(map(abs,pair)) == 2}
    assert set(held) == expected and all(v > center for v in held.values())
    assert all(v > max(3*center, arb("0.00050939")) for v in held.values())
    assert all(held[0,-2] < v for k,v in held.items() if k != (0,-2))
    for pair, cost in {**costs,**held}.items():
        print("编号偏移与最小须区分轨道的概率", pair, cost.str(22))



if __name__ == "__main__":
    verify()
