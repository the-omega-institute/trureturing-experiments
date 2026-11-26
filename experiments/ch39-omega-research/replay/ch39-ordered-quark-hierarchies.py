"""Two three-step ordered Gaussian-density minima with up and down quark hierarchies."""
from pathlib import Path
import runpy
from flint import arb,arb_poly,acb_mat

phase=runpy.run_path(str(Path(__file__).with_name("ch39-symmetric-source-transition.py")))
fit,finite=phase["fit"],phase["finite"]


def laws(budget):
    result=[]
    for orbit in range(4):
        det,k,(a,c,d)=phase["data"](orbit)
        a=1+(a-1)*budget/phase["B"]
        x=arb_poly([0,budget/3])
        y=budget-2*x
        p=a+c*(x*x+2*x*y)+d*x*x*y
        n=x*y*(c+d*x)
        result.append((det,k,p,n))
    return result


def select(budget,t0):
    candidates=laws(budget)
    det,k,p,n=candidates[3]
    tau=n(t0)/p(t0)
    def action(det,p,t):
        return (det/candidates[3][0]*p(t)).log()-tau*(t*t*(3-2*t)).log()
    expected=action(det,p,t0)
    gap=arb(10)
    for orbit,(base,_,poly,numerator) in enumerate(candidates):
        roots=phase["isolate"](numerator-tau*poly)
        assert all(0<t<1 for t in roots)
        for t in [arb(1)]+roots:
            if orbit==3 and t.contains(t0):
                continue
            diff=action(base,poly,t)-expected
            assert diff>0
            gap=min(gap,diff.lower())
    y=budget*(1-2*t0/3)
    x=budget*t0/3
    matrices=[]
    for pos in range(3):
        kk=acb_mat(k)
        for i in range(3):
            finite["add_edge"](kk,fit["inner"][i],fit["outer"][(i+1)%3],y if i==pos else x)
        degree,margin=finite["flux"]["degrees"](kk)
        assert degree==(3,[3]) and margin>0
        t=-fit["block"](kk,fit["ii"],fit["oo"])
        matrices.append(t*t*t)
    print("total stiffness, density exponent, isolated energy gap, weak and strong edges",
          *(v.str(10) for v in (budget,tau,gap,x,y)))
    return matrices


if __name__=="__main__":
    u=select(arb("4.160684"),arb(".0765543"))
    d=select(arb("24.600273"),arb(".0780958"))
    observed=[]
    for m in d:
        vals=finite["observables"](u[0],m)
        goals=finite["common"]["up"]+list(reversed(finite["common"]["down"]))
        assert all(abs(v/r-1)<arb(".0001") for v,r in zip(vals[:4],goals))
        observed.append([v.str(8) for v in vals[4:]])
    print("s12,s23,s13,delta,J at three relative orientations",observed)
