"""全部对角轻谱探针与一条相干辐边的前两个 RH 时间读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


root = Path(__file__).parent
model = run_path(str(root/"ch39-two-boundary-site-mixing.py"))
weights = run_path(str(root/"ch39-arithmetic-weight-rise.py"))
shell = model["shell"]


def profiles(numbers):
    sites = shell["geo"]["sites"]
    pole, tip = sites.index((17,17,17)), sites.index((13,23,13))
    labels = [i for i,(x,y,z) in enumerate(sites) if (x,y,z)==min((x,y,z),(y,z,x),(z,x,y))]
    rows, bonds, heights = [], [], []
    for n in numbers:
        h = acb.zeta_zero(n).imag
        neutral = shell["laplacian"](0,h)
        center = sum(neutral[k,k].real for k in range(19))/18
        heights.append(h)
        diagonal = []
        for q in (1,3,2,4,6):
            f = model["edge"]["mode_columns"](model["edge"]["filters"]["spectra"](q,h,center))
            diagonal.extend(sum(abs(f[s,k])**2 for k in range(3)) for s in labels)
            if q == 6:
                baseline = sum(abs(f[s,k])**2 for s in (pole,tip) for k in range(3))
                cross = 2*sum((f[tip,k]*f[pole,k].conjugate()).real for k in range(3))
                bonds.append((baseline,cross))
        rows.append(diagonal)
    return list(zip(*rows)),bonds,[(heights[0]/h)**2 for h in heights]


def normalize(row):
    assert all(v > 0 for v in row)
    total = sum(row)
    return [v/total for v in row]


def first_scale(w, nodes, first):
    lo,hi = arb(0),arb(1)
    for _ in range(110):
        mid = (lo+hi)/2
        a = sum(p*mid*u/(1+mid*u) for p,u in zip(w,nodes))
        if a < first:
            lo = mid
        elif a > first:
            hi = mid
        else:
            break
    return lo.union(hi)


def verify():
    ctx.prec,ctx.cap = 1536,3
    times = weights["joint"]["clock"]["times"]()
    first, second = 1-times[1], 1-times[1]-(1-times[2])/4
    ctx.prec = 320
    diagonals,bonds,u = profiles(weights["numbers"])
    high = arb(".004898")
    for row in diagonals:
        w = normalize(row)
        assert sum(p*high*z/(1+high*z) for p,z in zip(w,u)) > first
        b = sum(p*3*(high*z)**2/((1+high*z)*(1+4*high*z)) for p,z in zip(w,u))
        assert 4*(second-b) > arb("1e-5")
    def experiment(visibility):
        w = normalize([a-visibility*b for a,b in bonds])
        s = first_scale(w,u,first)
        diff = 4*(sum(p*3*(s*z)**2/((1+s*z)*(1+4*s*z)) for p,z in zip(w,u))-second)
        return diff,s
    lo,hi = arb(".99901809999"),arb(".99901810002")
    assert experiment(lo)[0] < 0 < experiment(hi)[0]
    w = normalize([a-lo.union(hi)*b for a,b in bonds])
    # 用两个有理端点同时包围整段能命中首项的尺度。
    low_s, high_s = arb(".00563198525"), arb(".00563198535")
    def first_at(s):
        return sum(p*s*z/(1+s*z) for p,z in zip(w,u))
    assert first_at(low_s)<first<first_at(high_s)
    s = low_s.union(high_s)
    third = sum(p/(1+9*s*z) for p,z in zip(w,u))
    gap = times[3]-third
    assert arb("4e-7") < gap < arb("5e-7")
    print("荷六最近辐边：命中前两项的可见度区间",lo.str(12),hi.str(12),"第三项亏额",gap.str(6))


if __name__ == "__main__":
    verify()
