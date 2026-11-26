"""同一批物理候选的正权算术时间读数与两个新增实际零点。"""
from itertools import combinations
from math import comb
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


root = Path(__file__).parent
clock = run_path(str(root / "ch39-arithmetic-height-bandwidth.py"))
physics = run_path(str(root / "ch39-chiral-joint-vacuum-density.py"))


def positive_witness(times, numbers, basis, y):
    degree = len(basis)-1
    ref = acb.zeta_zero(min(numbers)).imag
    nodes = {n:(ref/acb.zeta_zero(n).imag)**2 / arb(y) for n in numbers}
    target = clock["cleared_moments"](times,degree)
    # 给所有其余基点严格正的有理权，再对方阵求解。
    fixed = {n:arb("1e-7") for n in numbers if n not in basis}
    rhs = arb_mat([[target[j]-sum(w*nodes[n]**j for n,w in fixed.items())] for j in range(degree+1)])
    weights = arb_mat([[nodes[n]**j for n in basis] for j in range(degree+1)]).solve(rhs)
    tilted = fixed | {n:weights[i,0] for i,n in enumerate(basis)}
    unscaled = {}
    for n,w in tilted.items():
        for j in range(1,degree+1):
            w *= 1+j*j*nodes[n]/200
        unscaled[n] = w
    assert min(unscaled.values()) > 0
    assert sum(unscaled.values()).contains(1)
    for j in range(1,degree+1):
        read = sum(w/(1+j*j*nodes[n]/200) for n,w in unscaled.items())
        assert (read-times[j]).contains(0)
    print("读数、基点数、最小原概率",degree,len(nodes),min(unscaled.values()).str(10))


def exclude_next(times, numbers, degree):
    # 每个支撑点上的正多项式都必须给非负的矩泛函。
    heights = [acb.zeta_zero(n).imag for n in sorted(numbers)]
    sites = sorted(((heights[0]/h)**2 for h in heights), key=float)
    moments = clock["cleared_moments"](times,degree)
    facets = []
    for indices in combinations(range(len(sites)),degree):
        signs = {(-1)**sum(j>i for j in indices) for i in range(len(sites)) if i not in indices}
        if len(signs) != 1:
            continue
        p = [arb(signs.pop())]
        for j in indices:
            p = [v-sites[j]*w for v,w in zip([0,*p],[*p,0])]
        facets.append([r*m for r,m in zip(p,moments)])
    mean = moments[1]/moments[0]
    left, right = (sites[0]/mean).lower(), (sites[-1]/mean).upper()
    def negative(coeff, mid, radius):
        # 在中点作有限 Taylor 展开，避免窄窗上的符号抵消。
        deriv = [sum(coeff[j]*comb(j,k)*mid**(j-k) for j in range(k,degree+1))
                 for k in range(degree+1)]
        return deriv[0]+sum(abs(v)*radius**k for k,v in enumerate(deriv[1:],1)) < 0
    def cover(lo, hi, depth):
        mid, radius = (lo+hi)/2, (hi-lo)/2
        if any(negative(p,mid,radius) for p in facets):
            return 1
        assert depth < 35
        return cover(lo,mid,depth+1) + cover(mid,hi,depth+1)
    print("排除下一读数、支撑点数、尺度分段数",degree,len(sites),cover(left,right,0))


def verify():
    ctx.prec,ctx.cap = 1536,3
    times = clock["times"]()
    seven = [11661,11934,12207,22027,28325,48067,58214]
    positive_witness(times,seven,(11661,11934,22027,58214),"1.0143")
    moments = clock["cleared_moments"](times,4)
    c = arb(7)/10
    a,b,d = [moments[j+2]-2*c*moments[j+1]+c*c*moments[j] for j in range(3)]
    k = a*d/b**2
    span = k.sqrt()+(k-1).sqrt()
    heights = [acb.zeta_zero(i).imag for i in seven]
    assert max(heights)/min(heights) < arb("4.102") < arb("4.18") < span
    print("前四项所需高度比下界",span.str(13))
    positive_witness(times,seven+[8359],(8359,11661,11934,22027,58214),".60878")
    positive_witness(times,seven+[95478],(11661,11934,12207,22027,48067,95478),".95425")
    positive_witness(times,seven+[8359,95478],(8359,11661,12207,22027,28325,48067,95478),".54058")
    exclude_next(times,seven+[8359],5)
    exclude_next(times,seven+[95478],6)
    exclude_next(times,seven+[8359,95478],7)
    physics["actual_examples"]((8359,95478),tolerance=".0022",min_gap="1e-14",direct=True)


if __name__ == "__main__":
    verify()
