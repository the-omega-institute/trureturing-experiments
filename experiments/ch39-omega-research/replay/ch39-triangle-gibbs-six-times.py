"""原三窗、六个固定算术时间与单位荷三角作用的 Gibbs 逆温度。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
action = run_path(str(root/"ch39-mincut-triangle-action.py"))
angles = run_path(str(root/"ch39-mincut-three-window-angle.py"))
odds = run_path(str(root/"ch39-six-time-prior-odds.py"))
physical = angles["sign"]["physical"]


def harmonic_data():
    counts, coeff = action["triangles"](shell)
    assert sum(counts.values()) == 529
    base = sum(m-coeff[k].real for k,m in counts.items())
    return counts, coeff, base


def action_range(data):
    counts, coeff, base = data
    # 旧三窗的首个真实零点高度超过 8000；每条软化辐边最多引入 2δ/面。
    delta_cap = arb(23)/25*arb(2)**-20/8000**2
    error = sum(2*j*n for (_,_,j),n in counts.items())*delta_cap
    harmonics = {}
    for m in ((1,0),(0,1),(1,-1)):
        z = sum(c for (i,k,_),c in coeff.items() if (i,k)==m)
        harmonics[m] = z.real
        error += abs(z.imag)
    assert error < arb("1e-12")
    def excess(x,y):
        return sum(harmonics[m]*(1-t.cos()) for m,t in zip(harmonics,(x,y,x-y),strict=True))+arb(0,error)
    lower, upper = arb(".1531295"), arb(".4263063")
    curve_lo, curve_hi = angles["sign"]["density"]["eigen_band"]()
    lo, hi = arb("18.94088212"), arb("18.94414547")
    assert lo < curve_lo < curve_hi < hi
    six_cut = (15-arb(189).sqrt())/2+arb("3e-8")
    six = angles["cell_count"](6, six_cut)
    low, high = angles["cell_count"](2,lo), angles["cell_count"](2,hi)
    half = arb.pi()/16
    pending = [(1,0,0)]
    cells, max_n = 0,1
    while pending:
        n,i,j = pending.pop()
        cells += 1
        max_n = max(max_n,n)
        x,y = [((arb(2*k+1,1)/n)-1)*half for k in (i,j)]
        E = excess(x,y)
        if lower < E < upper:
            continue
        if (x*x+y*y).upper() < arb(".10124")**2 or abs(x-y) > arb(".215915"):
            continue
        if low(n,i,j)[0] >= 16 or high(n,i,j)[1] < 16 or six(n,i,j)[1] < 1:
            continue
        assert n < 2**26, (n,i,j,E)
        pending.extend((2*n,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
    print("基准作用、相对作用上下界、全域检验格数与最大二分",base.str(14),lower,upper,cells,max_n,flush=True)
    return lower,upper


def record(number, data, *, bits=100):
    counts, coeff, base = data
    ctx.prec = bits
    height = acb.zeta_zero(number).imag
    clocks = shell["geo"]["phases"](height)
    u,v = clocks[724],clocks[804]
    delta = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    energy = sum((1-delta)**j*(n-(coeff[i,k,j]*u**i*v**k).real)
                 for (i,k,j),n in counts.items())-base
    return height,energy


def temperatures(data, lower, upper):
    low_mass, _, tail_mass = odds["radial_certificate"](root)
    H = arb(13)*10**12
    tail_prior = ((H/(2*arb.pi())).log()+1)/(2*arb.pi()*H)+arb(".28")*(2*H.log()+arb(".5"))/H**2
    assert tail_prior < arb("3.59422e-13")
    for sign, number, floor in ((1,117930158,"97668.21"),(-1,71363114,"12079.05")):
        _, alpha, curve = physical["observables"](number)
        assert arb(".117") < alpha < arb(".119") and arb(".998") < curve < arb("1.002")
        height, energy = record(number,data)
        assert height < 10**8
        gap = energy-lower if sign==1 else upper-energy
        bound = (tail_mass/(low_mass*height**2*tail_prior)).log()/gap
        assert gap > 0 and bound > arb(floor)
        print("逆温度符号、低段编号、高度、作用增量、必要绝对值下界",
              sign,number,height.str(14),energy.str(13),bound.str(12),flush=True)


def finite_exclusion(data):
    pinned = run_path(str(root/"ch39-pinned-cp-zero-support.py"))
    ids = pinned["indices"]
    assert len(ids)==len(set(ids))==23
    # 已有名单的低、中、高编号块互相远离，使用有序编号与一次边界高度即可分类。
    tail = (143921293989490, *pinned["high_indices"])
    union = [*tail,142269938072,142270888235]
    limits = {i: record(i,data,bits=72) for i in union}
    idx = 13219432687559
    assert max(set(ids)-set(tail)) == idx
    assert record(idx,data,bits=72)[0] < 13*10**12
    assert all(i in ids and limits[i][0] > 13*10**12 for i in tail)
    left, right = (limits[i][1] for i in union[-2:])
    assert left < arb(".15815") and right > arb(".41608")
    assert all(arb(".16274") < limits[i][1] < arb(".40693") for i in tail)
    weight = sum(1/limits[i][0]**2 for i in tail)
    cap = weight*max(limits[i][0]**2 for i in union[-2:])
    assert cap < arb(".00000836")
    print("二十三点所有实逆温度的远尾概率上界",cap.str(12),flush=True)


if __name__ == "__main__":
    ctx.prec=160
    data=harmonic_data()
    extrema=action_range(data)
    temperatures(data,*extrema)
    finite_exclusion(data)
