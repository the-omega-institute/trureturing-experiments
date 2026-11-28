"""两个编号共用逐线十二个格点份额时，规范交汇与引力目标的相容性。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
SHELL = (2_306_888_100_085,2_306_888_100_086)
HEAVY = (2_345_666,2_345_668)
POLICIES = ((-28961,83986,34488,-100000),(-33046,95772,34502,-100000),
            (-31779,92072,34516,-100000),(-31517,91343,34502,-100000),
            (-31419,91026,34512,-100000))


def dot(u,v):
    return sum(a*b for a,b in zip(u,v,strict=True))


def read_data():
    def load(name):
        return run_path(str(ROOT/f"ch39-{name}.py"))
    shell = load("native-rh-ray-pvm")["shell_weights"]
    single, lib, backgrounds = load("two-zero-cutoff-gravity")["gravity_fields"](
        ROOT,HEAVY,tight_neutral=True)
    model,_,_,initial,_ = lib["chiral"]["source"]["mean"]["inputs"]()
    gauge = load("signed-gravity-unification")
    fields = gauge["gauge_weights"](model)
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    logscale = (arb(10**16)/model["mz"]).log()
    baseline = [a-b*logscale for a,b in zip(initial,beta,strict=True)]
    rows = []
    for code, background in zip(SHELL,backgrounds,strict=True):
        sector = {}
        for charge,levels,direct,_ in shell(ROOT,number=code):
            if charge in (-4,2):
                pivot = 0 if charge == -4 else 2
                sector[charge] = list(zip([z/levels[-1] for z in levels[:-1]]+[arb(1)],
                                          [w.real for w in direct[pivot]],strict=True))
        table = []
        for (s,w),(t,v) in product(sector[-4],sector[2]):
            left,right = s*t,(s+t)/2
            cell = []
            for k in range(12):
                x = arb(7)/1000*((11-k)*left+k*right)/11
                _,kg,wg,d1,d2,d3 = gauge["threshold_column"](
                    background,single["moments"],lib["expint"],fields,x)
                cell.append([w*v*a for a in (kg,wg,d1-d3,d2-d3,d3)])
            table.append(cell)
        rows.append(table)
    return rows,baseline,beta,model


def verify():
    rows,baseline,beta,model = read_data()
    gaps = [baseline[j]-baseline[2] for j in range(2)]*2
    slope = [beta[j]-beta[2] for j in range(2)]
    y = [arb(v)/10**9 for v in (-31975768,-1000000000,-667013,1938919,
                                  21636055,663870882,800655,-2324641)]
    lower = dot(y,[-1,-arb("1e-122"),*gaps[:2]]*2)
    selectors = [[arb(v)/10**5 for v in policy] for policy in POLICIES]
    combinations = [list(gaps) for _ in selectors]
    colour = [baseline[2] for _ in selectors]
    mismatch_gains = []
    for a,b in zip(*rows,strict=True):
        scores, observations = [], []
        first_scores, second_scores = [], []
        for first,second in zip(a,b,strict=True):
            first_scores.append(dot(y[:4],first[:4]))
            second_scores.append(dot(y[4:],second[:4]))
            scores.append(first_scores[-1]+second_scores[-1])
            observations.append(first[2:4]+second[2:4])
        lower += min(s.lower() for s in scores)
        # 两组在此线完全独立选格点，相对共用格点所能减少的投影。
        gain = (min(s.upper() for s in scores)-min(s.lower() for s in first_scores)
                -min(s.lower() for s in second_scores))
        mismatch_gains.append(gain.upper())
        for j,selector in enumerate(selectors):
            # 以中点选择格点即可；可行性由被选格点的完整区间响应验证。
            ranks = [(10**9*dot(selector, obs).mid()+arb(1)/2).mid().floor().unique_fmpz()
                     for obs in observations]
            k = min(range(12), key=ranks.__getitem__)
            combinations[j] = [g+x for g,x in zip(combinations[j],observations[k],strict=True)]
            colour[j] += b[k][4]
    lo = (model["mz"]/10**16).log()
    hi = (arb("1.22089e19")/10**16).log()
    for pair in (y[2:4],y[6:8]):
        speed = dot(pair,slope)
        lower -= (speed*lo).max(speed*hi)
    floor = lower/sum(abs(z) for z in y[2:4]+y[6:8])
    assert floor > arb("1.212")
    assert sum(mismatch_gains) > 0
    drift = lower / sum(mismatch_gains)
    total_drift = lower / max(mismatch_gains)
    assert drift > arb(".0072656") and total_drift > arb(".0079021")
    print("若两个背景各自精确匹配，逐线总变差距离的最大值与和的下界",
          drift.str(12),total_drift.str(12))
    matrix = arb_mat([v+[arb(1)] for v in combinations]).transpose()
    assert not matrix.det().contains(0)
    coeff = matrix.solve(arb_mat([[0],[0],[0],[0],[1]])).entries()
    assert all(v > arb(".00518") for v in coeff)
    min_colour = dot(coeff,colour)+beta[2]*(arb(10**16)/model["mz"]).log()
    assert min_colour < 0
    print("两组规范等式的五个布局混合份额",*[v.str(10) for v in coeff])
    print("该见证第二组在电弱标尺的色支，及引力精确时的规范最大残差下界",
          min_colour.str(12),floor.str(12))


if __name__ == "__main__":
    verify()
