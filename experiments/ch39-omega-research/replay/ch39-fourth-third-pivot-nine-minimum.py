"""用 1 TeV 的有限截止轻场核排除最后一种八层支撑。"""
from pathlib import Path
from runpy import run_path
from flint import arb


def light_bands(lo,cut,fine,p0,p1,lib):
    p=(arb(1000)/arb("91.1876"))**2
    rho=(arb("246.22")/arb("91.1876"))**2/2
    w=arb(".63811")**2*rho/2
    z=w+3*arb(".467453")**2*rho/10
    t=arb(".8616")**2*rho
    axial,boson=lib("neutral-axial-top"),lib("neutral-bosonic")
    pi=arb.pi()
    infinity=(-(2*axial["J"](w/p)+boson["j_neutral"](z/p))/(4*pi)
              +boson["scalar_difference"](w/p,z/p)/(8*pi)
              +3*axial["J"](t/p)/(8*pi)-axial["L"](t/p)/(24*pi))
    nangle=192
    endpoints=[arb(j)/nangle for j in range(nangle+1)]
    angles=[]
    for a,b in zip(endpoints,endpoints[1:]):
        ha,hb=a*(1-a),b*(1-b)
        angles.append((a,b,min(ha,hb),max(ha,hb)))
    def moments(x):
        def e1(y):
            return (y/x).expint(1)
        def e2(y):
            return arb(1) if y.is_zero() else (-y/x).exp()-y/x*e1(y)
        def block(m,h):
            return e1(m)-e1(m+p*h)
        def scalar(m,h):
            return x/p*(1-e2(p*h)-e2(m)+e2(m+p*h))
        sums=[arb(0) for _ in range(6)]
        def add(j,left,right):
            sums[j]+=left.lower().union(right.upper())
        for a,b,l,r in angles:
            for j,m in ((0,w),(2,t)):
                add(j,m/p*block(m,l),m/p*block(m,r))
            top=(p/z).log1p() if a.is_zero() else block(a*z,r)
            add(1,z/p*block(b*z,l),z/p*top)
            bottom=6*l*(e1(p*r)-e1(t+p*r))
            top=6*r*((t/(p*r)).log1p() if l.is_zero() else e1(p*l)-e1(t+p*l))
            add(3,bottom,top)
            add(4,scalar(a*z,l),scalar(b*z,r))
            add(5,scalar(w,l),scalar(w,r))
        return [v/nangle for v in sums]
    # 六个正角矩逐个随 x 增加，可以在相邻网格端点之间包围整个连续段。
    small,large=lo/p0,arb("2e5")
    axis=[small*((large/small).log()*j/320).exp() for j in range(321)]
    values=[moments(x) for x in axis]
    intervals=[]
    for left,right in zip(values,values[1:]):
        a,b,c,d,e,f=[u.lower().union(v.upper()) for u,v in zip(left,right,strict=True)]
        intervals.append((-(2*a+b)/4+(e-f)/8+3*c/8-d/24)/pi)
    kappa=arb("7.468e-36")
    span=(cut/lo).log()
    edges=[lo*(span*j/fine).exp() for j in range(fine+1)]
    bands=[]
    idx=0
    for l,r in zip(edges,edges[1:]):
        inverse=arb(0,(kappa/l).upper())
        if l/p0>axis[-1]:
            bands.append(inverse)
            continue
        while idx+1<len(intervals) and axis[idx+1]<l/p0:
            idx+=1
        enclosure=intervals[idx]
        j=idx+1
        while j<len(intervals) and not axis[j]>r/p0:
            enclosure=enclosure.union(intervals[j])
            j+=1
        correction=infinity-enclosure
        if not r/p0<axis[-1]:
            correction=correction.union(inverse)
        bands.append(correction)
    # 每个高端轻场正矩都由其无截止值控制；标量差分成两个正熵矩。
    high=lib("full-neutral-running")["weak_bound"]({"p0":p1},[p1],[arb(0)],axial,boson)
    return infinity,high,bands


def verify(root):
    core=run_path(str(root/"ch39-fourth-third-pivot-six-low-path-obstruction.py"))
    support=(1,3,10,15,16,17,18,19)
    probes=[(32000000,"y"),(100000000,"y"),(320000000,"y"),(560000000,"y"),
            (1000000000,"y"),(3160000000,"y"),(10000000000,"y"),
            (180000,"col"),(100000,"col"),(316000,"col"),(56000,"col"),(42000,"col"),
            (10000,"y"),(5600,"y"),(1000,"y")]
    shifts=[(p,ch,"-.22" if ch=="y" else ".58","0") for p,ch in probes]
    shifts.extend((p,"y","-.38",".5") for p in (320000000,560000000,1000000000,3160000000))
    _,_,constants,high_error,bounds=core["integral_table"](
        root,probes,shifts,resolution=(2048,32768),low_rows=support,precise_light=(1000,light_bands))
    limit=bounds.shape[1]//2-1
    core["integer_cover"](root,constants,high_error,bounds,8,[[(i,0,limit) for i in support]])


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
