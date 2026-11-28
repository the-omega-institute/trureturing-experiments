"""三条固定原线的共同逆耦合全范围、能标灵敏度及源耦合迁移。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
K = "4.12988"
TOL = "8e-7"


def source_coupling():
    """在源重谱、已冻结的布局及其 10¹² GeV 交点读出第三支。"""
    source = 2345698
    repair = run_path(str(ROOT/"ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    kernel,moments,_,_,fields,bgs,_,_,total,*_ = load("two-share-global-forward")["context"]((source,),20)
    model,_,_,initial,_ = load("chiral-minimal-neutral-gravity")["chiral"]["source"]["mean"]["inputs"]()
    coeff = kernel["series"](bgs[source],fields)
    third = arb_mat([[coeff[4,j] for j in range(coeff.ncols())]])
    n = kernel["ORDER"]+1
    tail = 2*arb(4)**n*total*sum(abs(a[2]) for a in fields)/(4*arb.pi()*factorial(n))
    return initial[2]+(third*moments)[0,0]+arb(0,tail.upper())-model["beta"][2]/(2*arb.pi())*(arb(10)**12/model["mz"]).log()


def verify():
    boundary = run_path(str(ROOT/"ch39-three-share-global-energy.py"))
    scale_top,scale_ew = boundary["verify"]()
    k,tol = arb(K),arb(TOL)
    point,mz = run_path(str(ROOT/"ch39-three-share-energy-interval.py"))["setup"](log_compensation=k)
    def read(s):
        return point(s,derivatives=True,common_readout=True)
    def bound(s):
        result = read(s)
        if result is None:
            return False
        _,_,dl,dc,remainder = result
        return dc>0 and dl<0 and abs(remainder)<tol*dc
    count = boundary["cover"](((scale_top.lower(),scale_ew.upper()),),bound)
    def centered(s):
        energy,common,dl,dc,_ = read(s)
        center = s.mid()
        e0,c0,*_ = read(center)
        return (energy.intersection(e0*((s-center)*dl).exp()),
                common.intersection(c0+(s-center)*dc))
    _,cmin = centered(scale_top)
    _,cmax = centered(scale_ew)
    upper = k+tol
    lower = k-tol
    affine_error = ((cmax-cmin)*tol).exp()-1
    one_percent = arb("1.01").log()/upper
    print("认证盒数、共同逆耦合的两个端点、全域指数式反演相对误差上界",
          count,cmin.str(12),cmax.str(12),affine_error.upper().str(10),flush=True)
    assert arb("35.574882")<cmin<arb("35.574884") and arb("43.118879")<cmax<arb("43.118881")
    assert one_percent>arb(".002409") and affine_error<arb("0.0000061")
    print("百分之一能标误差允许的共同逆耦合绝对误差下界",one_percent.lower().str(12),flush=True)
    cref = source_coupling()
    assert cmin<cref<cmax and arb("37.4959732")<cref<arb("37.4959733")
    scale = boundary["isolate"](scale_top,scale_ew,lambda s:read(s)[1]-cref)
    energy,_ = centered(scale)
    prediction = mz*(k*(cmax-cref)).exp()
    print("移用源端共同逆耦合时的耦合、标尺比、精确能标和近似反演（GeV）",
          *[v.str(13) for v in (cref,scale,energy,prediction)],flush=True)
    assert arb("1.10940e12")<energy<arb("1.10942e12") and abs(energy/prediction-1)<arb(".0000061")
    high = arb(10)**13
    crossing = boundary["isolate"](scale_top,scale_ew,lambda s:high-read(s)[0])
    _,c13 = centered(crossing)
    assert arb("36.963570")<c13<arb("36.963572")
    print("达到 10^13 GeV 的共同逆耦合阈值",c13.str(13),flush=True)
    print("负对数能标导数的统一上下界",lower.str(10),upper.str(10),flush=True)


if __name__=="__main__":
    verify()
