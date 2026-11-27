"""在最少八个活动站内，以两站门连续消去七壳弱差所需的最大平方弦长。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def components(sites, coordinates, maximum):
    groups = [{i} for i in sites]
    for a, b in combinations(sites, 2):
        squared = sum((x-y)**2 for x,y in zip(coordinates[a],coordinates[b],strict=True))
        if squared > maximum:
            continue
        left = next(group for group in groups if a in group)
        right = next(group for group in groups if b in group)
        if left != right:
            left.update(right)
            groups.remove(right)
    return groups


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    classification = lib("eight-optimal-classification")
    winners = classification["experiment"](root, fifth, columns, levels)
    graph = lib("three-ray-chiral-boundary")
    sites = graph["geo"]["sites"]
    trivial = lib("cyclic-two-characters-nine")["characters"](sites)[0][0]
    adj = classification["adj"]
    orbit = [[i for i in range(19) if abs(trivial[i,k]) > arb(".5")] for k in range(6)]
    floor, heavy, top = levels
    source = fifth*columns*acb_mat([[(heavy-floor).sqrt(), 0], [0, (top-floor).sqrt()]])
    means = adj(trivial)*source

    rows = []
    for j in range(7):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, 2_306_888_100_092+j)
        rows.append((adj(basis)*trivial).tolist()[0])
    ctx.prec=256
    caps = [((arb(q)/arb("1.22089e19"))-floor).sqrt()
            for q in (700,700,2000,2000,700,700,700)]
    charts = []
    for selected in combinations(range(7), 6):
        inverse = acb_mat([rows[j][:6] for j in selected]).inv()
        fixed = -inverse*acb_mat([[rows[j][6]] for j in selected])*acb_mat([means.tolist()[6]])
        charts.append((inverse,fixed,[caps[j] for j in selected]))

    def component_gap(part, groups):
        covered = [k for k, group in enumerate(groups) if part.issuperset(group)]
        if not covered:
            return None
        local = acb_mat([source.tolist()[i] for i in sorted(part)])
        gram = adj(local)*local
        scale = [(arb(3)/len(groups[k])).sqrt() for k in covered]
        base = acb_mat([[sum(source[i,h] for i in groups[k])/arb(len(groups[k])).sqrt()
                         -v*means[k,h] for h in (0,1)] for k,v in zip(covered,scale,strict=True)])
        gaps = []
        for inverse, fixed, cap in charts:
            J = acb_mat([inverse.tolist()[k] for k in covered])
            center = acb_mat([fixed.tolist()[k] for k in covered])
            gaps.append(classification["exclusion_test"](base, gram, scale, J, center, cap))
        return max(gaps)

    # 中继站可在门序中参与，最终在八站以外逐点恢复恒等。
    # 在 151 档一次放开其余全部十一站即可覆盖任意中继子集。
    bounds = ((243, 0, ".02"), (215, 1, ".02"), (151, 11, ".15"))
    min_gaps = {}
    connected = []
    for support in sorted(winners):
        active = set(support)
        spare = sorted(set(range(19)) - active)
        groups = [tuple(i for i in support if i in o) for o in orbit]
        if len(components(support,sites,244)) == 1:
            connected.append(support)
        cached = {}
        for cutoff, count, tolerance in bounds:
            for extra in combinations(spare,count):
                gaps = []
                for block in components((*support,*extra),sites,cutoff):
                    part = frozenset(block & active)
                    if part not in cached:
                        cached[part] = component_gap(part,groups)
                    if cached[part] is not None:
                        gaps.append(cached[part])
                assert gaps and max(gaps)>arb(tolerance), (support,extra,cutoff,gaps)
                gap = max(gaps)
                if cutoff not in min_gaps or gap < min_gaps[cutoff]:
                    min_gaps[cutoff] = gap
    assert len(winners)==32 and len(connected)==4
    print("不设中继时在平方弦长 244 连通的四种可行支撑", connected,flush=True)
    examples = ((216, (0,1,2,4,6,10,11,15), (16,)),
                (152, (1,3,4,5,8,9,10,12), (14,18)))
    for cutoff, support, extra in examples:
        assert support in winners and not set(support).intersection(extra)
        assert any(part.issuperset(support) for part in components((*support,*extra),sites,cutoff))
        print("可达平方弦长、最终八站、允许的中继站", cutoff, support, extra,flush=True)
    print("平方弦长上限与最小 Gram 违反量", {k:v.str(12) for k,v in min_gaps.items()},flush=True)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    source=run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root,*source["verify"](root,return_operator="spectral"))
