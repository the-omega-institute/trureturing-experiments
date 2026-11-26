"""十维几何刚度交点上的标准模型轻场两圈投影插入。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

vacuum = run_path(str(Path(__file__).with_name("ch39-ten-common-triangle-vacuum.py")))
flow = run_path(str(Path(__file__).with_name("ch39-two-loop-shell-separation.py")))
soft = vacuum["soft"]
crossing = vacuum["phase"]["crossing"]
sheet = crossing["smooth"]["sheet"]


def insertion_family(light,motion,r):
    shifts = arb(0,(arb(2).sqrt()*motion).upper())
    m = crossing["shell_spectra"]((0,5))
    m0 = [v+shifts for v in m[0]]
    m5 = [v+shifts for v in m[5]]
    mz,em,weak=arb("91.1876"),arb("127.930"),arb(".23129")
    initial = (3*(1-weak)*em/5,em*weak,1/arb(".1180"))
    pi=arb.pi()
    beta=flow["BETA"]
    p=(arb(115)/218,-arb(333)/218,arb(1))
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    assert dot(p,beta).contains(0) and sum(p).contains(0)
    gauge=[dot(p,col) for col in zip(*flow["B"])]
    yuk=[dot(p,col) for col in zip(*flow["D"])]
    assert all(v.contains(t) or v.overlaps(t) for v,t in zip(yuk,(arb(66)/109,-arb(3)/109,arb(3)/109)))
    adjoint = (1,3,-2)
    null = (5,-3,-2)
    common = (0,arb(2)/5,arb(3)/5)
    basis = (p,null,common)
    assert dot(null,adjoint)==0 and dot(common,adjoint).contains(0)
    # Two-loop Dirac coefficients for the three SM components of wedge^2 5.
    q = flow["q"]
    casimir_adjoint = (0,2,3)
    projected_dirac = []
    for indices,casimirs,levels in (
        (("4/5","0","1/2"),("4/15","0","4/3"),m5),
        (("1/10","3/2","1"),("1/60","3/4","4/3"),[light,*m0]),
        (("3/5","0","0"),("3/5","0","0"),m5)):
        T = list(map(q,indices))
        rows = [[4*dot(vv,T)*q(v)+q("20/3")*group*vv[j]*T[j]
                 for j,(v,group) in enumerate(zip(casimirs,casimir_adjoint))] for vv in basis]
        projected_dirac.append((rows,levels))

    def towers(Q2):
        def total(row,vector):
            # At small momentum the leading normalized positive log bound is
            # both faster and narrower than the cancelling closed formula.
            if Q2 < arb("1e-9"):
                return arb(0).union(sum(Q2/(8*v) for v in row).upper())
            return sum(sheet(v.upper(),Q2,vector).union(
                       sheet(v.lower(),Q2,vector)) for v in row)
        if Q2 < light.lower()*arb("1e-12"):
            first=arb(0).union((Q2/(8*light.lower())).upper())
        else:
            first=sheet(light.lower(),Q2,False).union(sheet(light.upper(),Q2,False))
        return total(m[0],True),total(m[5],True),total(m0,False)+first,total(m5,False)
    def heavy(v0,v5,d0,d5):
        return [7*(n*v0+k*v5)/(4*pi)-2*(zw*d0+fw*d5)/(3*pi)
                for n,k,zw,fw in zip((0,2,3),(5,3,2),(arb(".1"),arb("1.5"),1),(arb("1.4"),0,arb(".5")))]
    def insertion(Llo,Lhi,steps,yuk_cap=None):
        L=Llo.union(Lhi)
        ref2=r*r*(-2*Llo).exp()
        ref_err=arb(0,(ref2/8*(19*(35/(4*pi)+1/pi)+1/(pi*light.lower()))).upper())
        area=[arb(0) for _ in range(3)]
        heavy_range=[arb(0) for _ in basis]
        heavy_cost=arb(0)
        def endpoints(i):
            t=arb(i)/steps
            return [towers(r*r*(-2*(1-t)*x).exp()) for x in (Lhi,Llo)]
        left=endpoints(0)
        for i in range(steps):
            right=endpoints(i+1)
            block=[a.union(b) for a,b in zip(left[0],right[1])]
            time=arb(i).union(arb(i+1))*L/steps
            values=[a-b*time/(2*pi)+v+ref_err for a,b,v in zip(initial,beta,heavy(*block))]
            assert all(v>0 for v in values)
            max_Q2 = r*r*(-2*(1-arb(i+1)/steps)*Llo).exp()
            for rows,levels in projected_dirac:
                active_max = sum((max_Q2/v).min(arb(1)) for v in levels)
                for j,row in enumerate(rows):
                    coefficient = sum(entry/value for entry,value in zip(row,values))
                    heavy_range[j] += arb(0).union(-active_max*coefficient)
                    if j==0:
                        heavy_cost += active_max*coefficient.max(arb(0))
            for j,v in enumerate(values):
                area[j] += 1/v
            left=right
        gauge_insert=-dot(gauge,area)*L/(steps*8*pi**2)
        # Each down-type Yukawa eigenvalue squared may be as large as 4*pi.
        # Up and charged-lepton traces have positive coefficients.
        worst_yuk=yuk[1]*12*pi*Lhi/(32*pi**3)
        h_end=heavy(*towers(r*r))
        old=dot(p,[a+v+ref_err for a,v in zip(initial,h_end)])
        light_residual = old+gauge_insert+worst_yuk
        all_dirac = light_residual-heavy_cost.upper()*Lhi/(steps*8*pi**2)
        # Bound each Yukawa trace when enclosing all three projections.
        if yuk_cap is None:
            yuk_cap = 12*pi
        boxes = []
        for i,v in enumerate(basis):
            br = [dot(v,col) for col in zip(*flow["B"])]
            dr = [dot(v,col) for col in zip(*flow["D"])]
            one = dot(v,[a+t+ref_err for a,t in zip(initial,h_end)])-dot(v,beta)*L/(2*pi)
            ybox = sum(arb(0).union(t*yuk_cap) for t in dr)*L/(32*pi**3)
            boxes.append(one+(heavy_range[i]-dot(br,area))*L/(steps*8*pi**2)+ybox)
        return gauge_insert,light_residual,all_dirac,boxes
    return insertion,mz,m[5],dot(p,adjoint)


def verify():
    ctx.prec=320
    h = acb.zeta_zero(11431).imag
    r = soft["geo"]["theta"](804)/soft["geo"]["theta"](12)
    K = (h*h+soft["point"]["study"]["offset"])/r**2
    op,_,vnorm,squared,w = vacuum["setup"](h)
    light,motion = vacuum["light"](K,op,vnorm,squared,w)
    insertion,mz,vectors,projection = insertion_family(light,motion,r)
    low,high = ((arb(mu)/mz).log() for mu in ("1e13","1e17"))
    # A finite cover of the matching-scale window avoids binding every
    # early- and late-momentum sample to opposite endpoints of four decades.
    results=[]
    for j in range(8):
        a=low+(high-low)*j/8
        b=low+(high-low)*(j+1)/8
        results.append(insertion(a,b,64))
    gauge_low=min(v.lower() for v,_,_,_ in results)
    gauge_high=max(v.upper() for v,_,_,_ in results)
    residual=min(v.lower() for _,v,_,_ in results)
    print("SM gauge-only two-loop insertion range; residual lower including Yukawa bound",
          gauge_low.str(10),gauge_high.str(10),residual.str(10))
    assert gauge_low > arb(".685") and gauge_high < arb(".845") and residual > arb(".654")
    all_dirac = min(v.lower() for _,_,v,_ in results)
    print("Residual lower with all allowed Dirac decoupling functions",all_dirac.str(10))
    assert all_dirac > arb(".4617") and all_dirac*218/333 > arb(".3022")

    return insertion,mz,r,vectors,projection


if __name__ == "__main__":
    verify()
