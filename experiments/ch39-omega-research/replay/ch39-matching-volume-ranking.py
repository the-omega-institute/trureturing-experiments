"""Compare leading matching-volume weights in one common high-energy Yukawa box."""
from pathlib import Path
from runpy import run_path
from flint import arb, acb, ctx

root = Path(__file__).parent
gauge = run_path(str(root/"ch39-floating-scale-winner.py"))
weights = run_path(str(root/"ch39-hardy-height-selection.py"))
q = gauge["q"]


def verify():
    ctx.prec = 144
    top,interior = 10875,11488
    full = q(".999").union(q("1.001"))
    low_y = [2*(c*q(".99").union(q("1.01"))).log() for c in gauge["yukawa_centers"]]
    top_lo,top_hi = q("26.859461"),q("26.859504")
    for time,i in ((top_lo,1),(top_hi,2)):
        inverse = [1/(c*full)**2 for c in gauge["centers"]]
        end,_,_ = gauge["evolve"](inverse,low_y,time,16384)
        difference = end[i]/gauge["ratios"][i]-gauge["target"](top,time)
        assert (difference<0 if i==1 else difference>0)
    def rate_check(step,x,rate):
        if step >= 23*64:
            assert rate[1]>0 and rate[2]>0
    gauge["running"]["integrate"]([c*full for c in gauge["centers"]],33*64,rate_check)

    common_logs = [q(v)+arb(0,"1e-8") for v in
                   "-25.429 -12.978 -1.537 -23.804 -17.827 -10.108 -25.635 -14.934 -9.268".split()]
    loop = gauge["loop"]
    ratios = gauge["ratios"]
    beta,b,d = [gauge["running"][k] for k in ("BETA","B","D")]
    traces = [sum(y.exp() for y in common_logs[j:j+3]) for j in (0,3,6)]

    def density(n,lo,hi,steps,require_feasible):
        time = lo.union(hi)
        value = gauge["target"](n,time)
        terminal = [value*r for r in ratios]
        log_jac = arb(0)
        def observe(step,g,tr):
            nonlocal log_jac
            div = 2*sum(b[i][i]*v*v for i,v in enumerate(g))/loop**2
            div += (9*(tr[0]+tr[1])+5*tr[2])/loop
            log_jac += step*div
        x,y,_ = gauge["evolve"](terminal,common_logs,-time,steps,diagonal_quarks=True,observe=observe)
        if require_feasible:
            assert all(q(".999")<1/(v.sqrt()*c)<q("1.001") for v,c in zip(x,gauge["centers"]))
            assert all(q(".99")<(v/2).exp()/c<q("1.01") for v,c in zip(y,gauge["yukawa_centers"]))
        rates = []
        for i,r in enumerate(ratios):
            correction = sum(bb/v-dd*t for bb,dd,v,t in zip(b[i],d[i],terminal,traces))
            rates.append(2*value-2*(beta[i]+correction/loop)/(r*loop))
        assert all(v>0 for v in rates)
        return 4*ratios[0]*ratios[1]*ratios[2]*log_jac.exp()*sum(rates)

    a,bound = q("26.81235"),q("26.81262")
    top_volume = (top_hi-top_lo)*density(top,top_lo,top_hi,16384,False)
    lower = (bound-a)*density(interior,a,bound,16384,True)
    density(top,q("26.85948205"),q("26.85948210"),131072,True)
    ratio = lower/top_volume
    w_top,w_other = [weights["weight"](acb.zeta_zero(n).imag) for n in (top,interior)]
    weighted = ratio*w_other/w_top
    assert w_other<w_top and ratio>q("6.25") and weighted>q("5")
    print("11488/10875 matching volume >",ratio.lower().str(10),"including spectral weight >",
          weighted.lower().str(10),"spectral weight ratio",(w_other/w_top).str(10))


if __name__=="__main__":
    verify()
