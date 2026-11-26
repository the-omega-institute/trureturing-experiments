"""有限 RH Gram 探针消去一项弱差时，最低截止能量的锐上限。"""
from pathlib import Path
from runpy import run_path
from flint import arb, ctx


def ball(lo, hi):
    return arb(str(lo)).union(arb(str(hi)))


def verify(root):
    groups = run_path(str(root/"ch39-rh-cross-spectral-cutoff.py"))["verify"]()
    rows = [row for family in groups.values() for row in family]
    ordered = []
    for i in (1,2,3):
        a,b=(i,()),(i,tuple(range(i)))
        pair,=[row for row in rows if set(row[2:])=={a,b}]
        ordered.append((f"顺序枢轴{i}",pair))
    budgets = [(key, max(groups[key],key=lambda row:row[0])) for key in ("原向量","同一非空消元集")]
    budgets.extend(ordered)
    budgets.append(("全体",max(rows,key=lambda row:row[0])))
    ctx.prec = 192
    mz,planck=arb("91.1876"),arb("1.22089e19")
    p0,p1=[(m/planck)**2 for m in (mz,arb("1.176397e14"))]
    rho=(arb("246.22")/mz)**2/2
    g1,g2,yt=map(arb,(".467453",".63811",".8616"))
    w,z,t=[rho*m for m in (g2**2/2,(g2**2+3*g1**2/5)/2,yt**2)]
    high=-run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"](arb(1),p0,p1)
    weak=run_path(str(root/"ch39-neutral-mz-cutoff.py"))["low_weak_column"]
    def response(energy):
        return -weak(w,z,t,p0,p1,x_over_p0=(energy/mz)**2)
    intervals = [(7,8),(8,9),(181,182),(510,514),(1550,1560),(2650,2675)]
    for (name,row),(lower,upper) in zip(budgets,intervals,strict=True):
        b=ball(*row[:2])
        target=high*b/(1+b)
        lo,hi=arb(lower),arb(upper)
        assert response(lo)<target<response(hi)
        while hi-lo>arb(".02"):
            mid=(lo+hi)/2
            value=response(mid)
            if value<target:
                lo=mid
            else:
                assert value>target
                hi=mid
        print(name,"截止能量上限 GeV",lo.str(13),hi.str(13),flush=True)
    fields=run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))["gravity_fields"]
    index=2_306_888_100_085
    _,_,spectra=fields(root,(index-1,index,index+1))
    for number,spectrum in zip((index-1,index,index+1),spectra,strict=True):
        smallest = min(z.lower() for _,_,z in spectrum[1:])
        largest = max(z.upper() for _,_,z in spectrum[1:])
        assert smallest>arb("1.6e-11") and largest<1
        error=arb("2.398e-21")+arb("7.468e-36")*(1+2*arb("3844.11"))/smallest
        assert error<arb("7e-21")
        print(number,"重谱固有截止下的最小能量、全体实交叉读数的绝对误差界",
              (planck*smallest.sqrt()).lower().str(12),error.upper().str(12))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
