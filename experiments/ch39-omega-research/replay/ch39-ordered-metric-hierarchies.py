"""Common Maxwell geometry, both mass hierarchies and three wrong CKM branches at zero 10."""
from pathlib import Path
import runpy
from flint import arb,acb_mat,arb_poly

phase=runpy.run_path(str(Path(__file__).with_name("ch39-symmetric-source-transition.py")))
metric=runpy.run_path(str(Path(__file__).with_name("ch39-full-metric-ckm.py")))
fit,finite=phase["fit"],phase["finite"]
models=[phase["data"](i) for i in range(4)]
params=[
    ("4.407162902 -.081641631 .000675049 -.000029672 .000001672 -.000000106",
     ".07703544371 .00003073400 -.00003244799 .00000047893 -.00000004632 .00000000279"),
    ("23.14251745 .42286971 .007062482 -.000042563 .000001843 -.000000035",
     ".07857148143 -.00011299327 -.00000675697 .00000005980 -.00000000109 .00000000002")]
N=512
step=(metric["qhi"]-metric["qlo"])/N


def candidates(q,sector):
    u=(q-arb(".89"))/arb(".035")
    cheb=[arb(1),u]
    for _ in range(4):
        cheb.append(2*u*cheb[-1]-cheb[-2])
    budget,t0=[sum(x*arb(c) for x,c in zip(cheb,poly.split())) for poly in params[sector]]
    det,k,(a,c,d)=models[3]
    a=1+budget*(a-1)/phase["B"]
    for other,_,(aa,cc,dd) in models[:3]:
        aa=1+budget*(aa-1)/phase["B"]
        scale=other/det
        low=scale*aa-a+min(arb(0),(scale*cc-c).lower())*budget**2/3
        low+=min(arb(0),(scale*dd-d).lower())*(budget/3)**3
        assert low>0
    x=arb_poly([0,budget/3])
    y=budget-2*x
    p=a+c*(x*x+2*x*y)+d*x*x*y
    n=x*y*(c+d*x)
    tau=n(t0)/p(t0)
    roots=phase["isolate"](n-tau*p)
    assert sum(r.contains(t0) for r in roots)==1
    def energy(t):
        return p(t).log()-tau*(t*t*(3-2*t)).log()
    for t in [arb(1)]+[r for r in roots if not r.contains(t0)]:
        assert energy(t)>energy(t0)
    D=acb_mat([[q,0,0],[0,1,0],[0,0,1]])
    mats=[]
    for pos in range(3):
        m=acb_mat(k)
        for i in range(3):
            val=budget*(1-2*t0/3 if i==pos else t0/3)
            finite["add_edge"](m,fit["inner"][i],fit["outer"][(i+1)%3],val)
        assert finite["flux"]["degrees"](m)[0]==(3,[3])
        T=-fit["block"](m,fit["ii"],fit["oo"])
        mats.append(D*metric["adj"](metric["F"])*(T*T*T)*metric["F"]*D)
    return mats


if __name__=="__main__":
    metric["metric_interval"]()
    groups=[finite["common"]["up"],list(reversed(finite["common"]["down"]))]
    for i in range(N):
        q=metric["qlo"]+(arb(i)+arb(1)/2)*step
        u=candidates(q,0)[0]
        for d in candidates(q,1):
            vals=finite["observables"](u,d)
            assert all(abs(a/b-1)<arb(".00001") for a,b in zip(vals[:4],groups[0]+groups[1]))
            assert vals[0]/vals[1]<arb(".052") and vals[2]/vals[3]<arb(".052")
            assert vals[6]>arb(".122")
    dilation=(1+step/(2*metric["qlo"]))**2
    assert dilation*arb("1.00001")<arb("1.00018")
    slope=arb(2).sqrt()/metric["qlo"]*(1+arb(".06"))/(1-arb(".06"))
    assert arb(".122")-slope*step/2>arb(".1218")
    print("all physical metric values, 4 mass-ratio relative errors < .00018, every |Vub| > .1218")
