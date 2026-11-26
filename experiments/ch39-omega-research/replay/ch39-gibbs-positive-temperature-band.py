"""六时间必需窄段中，实际三窗零点的低作用隔离。"""
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx

common = run_path(str(Path(__file__).with_name("ch39-gibbs-negative-temperature-band.py")))
old, angles, shell = [common[k] for k in ("old","angles","shell")]
lower = arb(".1531305")
target_sum, target_diff = arb(".14318178"), arb(".00030820")
sum_radius, diff_radius = arb(".000003"), arb(".0007")


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
        if any(abs(x+y-s*target_sum)<sum_radius and abs(x-y-s*target_diff)<diff_radius for s in (-1,1)):
            continue
        energy = sum(w*(1-t.cos()) for w,t in zip(weights,(x,y,x-y)))
        if energy-arb("1e-12") > lower:
            continue
        if abs(x-y)>arb(".215915") or low(n,i,j)[0]>=16 or high(n,i,j)[1]<16 or six(n,i,j)[1]<1:
            continue
        assert n<2**28,(n,i,j,x,y)
        pending.extend((2*n,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
    print("正端作用的和差相位方框、角格数",sum_radius.str(9),diff_radius.str(9),cells,flush=True)



def screen_height(intervals):
    # 和差模 2π 会包括两钟同在 π 附近的提升；其余区间再检验真实时钟直线上的曲率与作用。
    a,b = shell["geo"]["a"],shell["geo"]["b"]
    _, coeff, _ = old["harmonic_data"]()
    harmonics = [sum(c for (i,j,_),c in coeff.items() if (i,j)==m).real
                 for m in ((1,0),(0,1),(1,-1))]
    cutoff=arb("18.94414547")
    fixed, terms=angles["inertia"]["schur_coefficients"](2,cutoff)
    kept=[]
    for left,right in intervals:
        possible=[]
        stack=[(left,right,0)]
        while stack:
            l,h,k=stack.pop()
            I=l.union(h)
            X,Y=a*I,b*I
            if X.cos() < (arb.pi()/16).cos():
                continue
            actual=sum(w*(1-t.cos()) for w,t in zip(harmonics,(X,Y,X-Y),strict=True))
            if actual-arb("1e-12") > lower:
                continue
            u,v=acb(0,2*X).exp(),acb(0,2*Y).exp()
            # 必要条件为 λ₁₆ < cutoff，故负惯性数至少 16。
            max_neg=16-fixed
            for c,aa,bb,cc in terms:
                s=c+2*(aa*u+bb*v+cc*u.conjugate()*v).real
                max_neg+=int(not(s>0))
            if max_neg<16:
                continue
            if k==6:
                possible.append(I)
                continue
            mid=(l+h).mid()/2
            stack.extend(((l,mid,k+1),(mid,h,k+1)))
        if possible:
            kept.append((min(z.lower() for z in possible),max(z.upper() for z in possible)))
    return kept


def lookup():
    a,b = shell["geo"]["a"],shell["geo"]["b"]
    intervals = []
    for sign in (-1,1):
        intervals.extend(common["physical_times"](a+b,b-a,common["low_h"],common["high_h"],
                         sign*target_sum,-sign*target_diff,sum_radius,diff_radius))
    print("正端和差时钟的连续高度区间数",len(intervals),flush=True)
    intervals=screen_height(intervals)
    print("相位提升、曲率与三角作用联合筛选后的区间数",len(intervals),flush=True)
    direct=0
    count=0
    for left,right in intervals:
        # 原研究的 RH 输入下，临界线上 ζ 的区间不含零即可排除。
        if not acb(arb(1)/2,left.union(right)).zeta().contains(0):
            direct+=1
            continue
        nl=int(left.zeta_nzeros().unique_fmpz())
        nr=int(right.zeta_nzeros().unique_fmpz())
        assert nr>=nl
        for n in range(nl+1,nr+1):
            count+=1
            h,e=old["record"](n,old["harmonic_data"](),bits=96)
            assert e>lower,(n,e)
    print("直接区间排除数、区间内实际零点数",direct,count,flush=True)


if __name__=="__main__":
    ctx.prec=112
    isolate_angles()
    lookup()
