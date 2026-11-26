"""由六时间所必需的窄高度段限制负逆温度尾端。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, ctx


root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
old = run_path(str(root/"ch39-triangle-gibbs-six-times.py"))
angles = old["angles"]
shell = old["shell"]
low_h, high_h = 41339000000, 41399000000
threshold = arb(".426303")
radius = arb(".000030")
target_x, target_y = arb(".16070324"), arb("-.05199106")


def required_band():
    heights = run_path(str(root/"ch39-six-time-height-bands.py"))
    clock = run_path(str(root/"ch39-arithmetic-height-bandwidth.py"))
    ctx.prec = 384
    times = clock["times"]()[:7]
    times[6] += arb("1.360900010e-8")
    denom, moments = heights["radial"]["compact_moments"](times)
    for h in (low_h,high_h):
        point=arb(h)
        left = (point-arb("1e-6")).zeta_nzeros().unique_fmpz()
        right = (point+arb("1e-6")).zeta_nzeros().unique_fmpz()
        assert left is not None and left == right
    bounds = []
    for h, upper in ((low_h,False),(high_h,True)):
        cutoff = 1/(1+(arb(h)/heights["scale"])**2)
        nodes, weights = heights["cutoff_rule"](cutoff,moments,denom)
        heights["verify_indicator"](nodes,denom,upper)
        tail = sum(w for t,w in zip(nodes,weights) if t<cutoff)
        bounds.append(tail+weights[0] if upper else tail)
    assert bounds[0]-bounds[1] > arb(".10448811")
    print("必需窄段的概率下界",(bounds[0]-bounds[1]).str(12),flush=True)


def isolate_angles():
    _, coeff, _ = old["harmonic_data"]()
    modes = ((1,0),(0,1),(1,-1))
    weights = [sum(c for (i,j,_),c in coeff.items() if (i,j)==m).real for m in modes]
    half = arb.pi()/16
    six = angles["cell_count"](6,(15-arb(189).sqrt())/2+arb("3e-8"))
    low, high = (angles["cell_count"](2,arb(t)) for t in ("18.94088212","18.94414547"))
    pending = [(1,0,0)]
    cells = 0
    while pending:
        n,i,j = pending.pop()
        cells += 1
        x,y = [(arb(2*k+1,1)/n-1)*half for k in (i,j)]
        if any(abs(x-s*target_x)<radius and abs(y-s*target_y)<radius for s in (-1,1)):
            continue
        energy = sum(w*(1-t.cos()) for w,t in zip(weights,(x,y,x-y)))
        if energy+arb("1e-12") < threshold:
            continue
        if abs(x-y)>arb(".215915") or low(n,i,j)[0]>=16 or high(n,i,j)[1]<16 or six(n,i,j)[1]<1:
            continue
        assert n<2**25,(n,i,j,x,y)
        pending.extend((2*n,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
    print("负端作用的两相位方框、角格数",radius.str(9),cells,flush=True)


def integer_part(x, method):
    return int(getattr(x,method)().unique_fmpz())


def physical_times(a,b,lo,hi,x,y,r,second_radius=None):
    s = r if second_radius is None else second_radius
    # 第一只钟写成 ah=2πk+x+e，|e|≤r，再按 b/a 的有理逼近列尽所有 k。
    tau = 2*arb.pi()
    ratio = b/a
    first = integer_part((a*lo-x-r)/tau,"ceil")
    last = integer_part((a*hi-x+r)/tau,"floor")
    approx = Fraction(str(float(ratio))).limit_denominator(10**11)
    p,q = approx.numerator,approx.denominator
    mod = pow(p,-1,q)
    target = (y-ratio*x)/tau
    # 只在 k-first 的范围积累有理逼近误差。
    shifted_target = target-(ratio*first)
    width = (s+r*abs(ratio))/tau+(last-first)*abs(ratio-arb(p)/q)
    lower = integer_part(q*(shifted_target-width),"ceil")
    upper = integer_part(q*(shifted_target+width),"floor")
    assert last-first<q
    seen = set()
    for t in range(lower,upper+1):
        j = t*mod % q
        if j>last-first:
            continue
        k = first+j
        assert k not in seen
        seen.add(k)
        h = (tau*k+x)/a
        bcenter = (b*h-y)/tau
        ell = integer_part(bcenter+arb(1)/2,"floor")
        h2 = (tau*ell+y)/b
        starts = (h-r/a).max(h2-s/b).max(lo)
        ends = (h+r/a).min(h2+s/b).min(hi)
        if starts<ends:
            yield starts.lower(),ends.upper()
        else:
            assert starts>ends


def lookup():
    a,b = shell["geo"]["a"],shell["geo"]["b"]
    intervals=[]
    for s in (-1,1):
        intervals.extend(physical_times(a,b,low_h,high_h,s*target_x,s*target_y,radius))
    count=0
    for left,right in intervals:
        nl=int(left.zeta_nzeros().unique_fmpz())
        nr=int(right.zeta_nzeros().unique_fmpz())
        assert nr>=nl
        for n in range(nl+1,nr+1):
            count+=1
            print("角盒可能实际零点",n,flush=True)
            # 如有根，直接检验模型作用或窗口。
            h,e=old["record"](n,old["harmonic_data"](),bits=96)
            assert e < threshold
    print("覆盖连续高度区间、其中实际零点数",len(intervals),count,flush=True)


if __name__=="__main__":
    required_band()
    ctx.prec=128
    isolate_angles()
    lookup()
