"""由全部随字长不增的字概率到 7.2% 末圈门槛的总变差。"""

from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


POSITIVE_DIRECTION = (
    5219145,-5554600,1112549,-92518,168800,618357,976269,195636,461617,19903,
    401441,-326133,720244,-130310,722163,-32713,731246,-101364,585435,-17313,
    554735,-43801,627507,-53960,599083,-82953,649856,-107794,661083,-42540,
    633612,-52986,609480,-51411,562242,-36338,599604,-48110,622380,-37229,
    618632,-62591,622616,-52121,588906,-24786,583388,-25142,568387,-6125,
    565818,-19027,595516,-18560,586880,-15887,583124,-28345,553346,5891,
    531296,22653,547153,41444,530565,36771,541558,1273,537986,27528,
    523544,33690,517159,62346,454934,107838,470027,110301,479880,132767,
    445634,51839,457437,52273,410952,186471,445160,225748,350196,324112,
    161852,266043,295448,279188,224953,374471,211497,8181,55883,1037972,
    -524258,2959134,1483501,3603257,
)
DESTINATIONS = {3:7917321, 5:4216352, 80:52144, 82:4789970, 84:7376239,
                86:7685693, 88:5555808, 90:756525}


def moments(matrix, p):
    return (matrix*arb_mat([[v] for v in p])).entries()


def unit_scores(rows, direction, sign):
    norm = sum(v*v for v in direction).sqrt()
    slope = (arb(9)/116).sqrt()
    return [sign*rows[-1][j] - slope*sum(v*r[j] for v, r in zip(direction, rows[:-1]))/norm
            for j in range(107)]


def lower_bound(score):
    running = score[0]
    means = [running]
    for r in range(1, 107):
        running += 2*score[r]
        means.append(running/(2*r+1))
    support = max(t.upper() for t in means)
    oscillation = max(t.upper() for t in score) - min(t.lower() for t in score)
    assert support < 0
    return -support/oscillation


def direct(rows):
    dims = [1]+[2]*106
    lambdas = {1:arb(215211930)/10**12, 105:arb(325198473300)/10**12}
    lambdas[106] = 1-sum(lambdas.values())
    mono = [d*sum(v/(2*r+1) for r, v in lambdas.items() if r >= j) for j, d in enumerate(dims)]
    change = [arb(DESTINATIONS.get(j, 0))/10**12 for j in range(107)]
    change[0] = -sum(change[1:])
    crossing = [p+arb("1.000001")*d for p, d in zip(mono, change)]
    matrix = arb_mat(rows)
    raw = [moments(matrix, p) for p in (mono, crossing)]
    powers = []
    for m in raw:
        total = sum(z*z for z in m)
        assert total > 0
        powers.append([v*v/total for v in m])
    assert all((sum(p)-1).contains(0) and min(v/d for v,d in zip(p,dims)) > arb(".003167")
               for p in (mono,crossing))
    assert arb(".06362451") < powers[0][-1] < arb(".06362453")
    assert arb(".072") < powers[1][-1] < arb(".07200002") and raw[1][-1] < 0
    positive = lower_bound(unit_scores(rows, list(map(arb, POSITIVE_DIRECTION)), 1))
    # 使用未作严格可行放大的方向，保留更紧的负支对偶。
    trial = moments(matrix, [p+d for p,d in zip(mono,change)])
    negative = lower_bound(unit_scores(rows, trial[:-1], -1))
    upper = -change[0]*arb("1.000001")
    assert positive > arb(".00072492") and negative > arb(".00003835004803")
    assert upper < arb(".00003835009036")
    print("正支下界、负支下界、过线见证的移动量",
          *(v.str(18) for v in (positive,negative,upper)))
    print("单调与过线两分布的第105圈份额",*(w[-1].str(18) for w in powers))
    return powers


def physical(root, powers):
    load = lambda name: run_path(str(root/("ch39-"+name+".py")))
    winding = load("mincut-global-winding")
    shell = load("three-ray-chiral-boundary")
    local = load("mincut-reference-ckm-action")
    flavor = load("chiral-neutral-flavor-phase")
    tri = winding["terms"](shell)
    distance = -local["reference_angles"](flavor)[0]
    offset = shell["soft"]["point"]["study"]["offset"]
    assert acb.zeta_zero(1).imag > 14 and offset > 0
    cap = arb(23)/25*arb(2)**-20/14**2
    idx = range(1, 106)

    def modes(delta, weights):
        return [(q,[weights[q-1]*v for v in winding["coeff"](tri,q,delta)]) for q in idx]

    for weights, floor, ceiling in zip(powers, (334,354), (467,486)):
        all_heights = modes(arb(0).union(cap), weights)
        slope = []
        for k in range(16):
            t = arb(".00815")+arb(".00465")*arb(k).union(arb(k+1))/16
            slope.append(winding["derivatives"](all_heights,-t,distance-t)[0])
        assert min(d.lower() for d in slope) > floor and max(d.upper() for d in slope) < ceiling
        print("全部正零点高度和合格参考味角段上的 x 偏导",floor,ceiling)

    for number, ranges in ((11661, ((-28,-27),(-42,-41))),
                           (2_306_888_100_085, ((409,410),(429,430)))):
        h = acb.zeta_zero(number).imag
        phases = shell["geo"]["phases"](h)
        x,y = [phases[d].arg() for d in (724,804)]
        delta = arb(23)/25*arb(2)**-20/(h*h+offset)
        readings = []
        for w, (a,b) in zip(powers,ranges):
            dx = winding["derivatives"](modes(delta,w),x,y)[0]
            assert a < dx < b
            readings.append(dx.str(15))
        print("实际零点编号、两钟角、两个作用的 x 偏导",number,x.str(14),y.str(14),*readings)


def differences(root, rows):
    diff = [[b-a for a,b in zip(row,after)] for row,after in zip(rows,rows[1:])]
    # 中心字 U0 给零向量。对其他中心均匀段，以锥投影证实非零差分不能过线。
    centers = []
    for row in diff:
        running = row[0]
        col = [running]
        for j in range(1,107):
            running += 2*row[j]
            col.append(running/(2*j+1))
        centers.append(col)
    for sign,active in ((-1,tuple(range(1,106,2))), (1,(*range(2,63,2),105,106))):
        X = arb_mat([[row[i] for i in active] for row in centers])
        last = arb_mat([[sign*centers[-1][i]] for i in active])
        coeff = (X.transpose()*X).solve(last)
        proj = (X*coeff).entries()
        assert all(x>0 for x in coeff.entries()) and sum(x*x for x in proj) < arb(".061")
        residual = [arb(sign*int(i==104))-v for i,v in enumerate(proj)]
        for j in range(1,107):
            if j not in active:
                assert sum(h*r[j] for h,r in zip(residual,centers)) < 0
    # 已有严格越线分布；其与 U0 混合后差分逐项只乘同一正数。
    source = run_path(str(root/"ch39-word-mass-variation.py"))
    change = [arb(0)]*107
    for j,d in source["TV_SHIFTS"]["diff",1]:
        change[j] = arb(d)/10**12
    change[-1] = -sum(change[:-1])
    p = [arb(n)/213+arb("1.0000002")*v for n,v in zip([1]+[2]*106,change)]
    m = moments(arb_mat(diff),p)
    share = m[-1]**2/sum(v*v for v in m)
    assert min(p)>0 and share > arb(".072")
    print("掺入中心字时保持不变的差分末圈份额",share.str(17))


def verify(root):
    geometry = run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))
    ctx.prec = 256
    theta = geometry["unit"]((108469,109)).arg()
    rows = [[(q*j*theta).cos() for j in range(107)] for q in range(1,107)]
    powers = direct(rows[:-1])
    physical(root,powers)
    differences(root,rows)


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
