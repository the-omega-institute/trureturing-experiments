"""原第三枢轴匹配两个引力矩时，分别限定正权和负权截止的下界。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number=2_306_888_100_089
    _,ordered=lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=number-3)[2]
    _,_,pivots,_=next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root,number) if row[0]==6)
    c=[v.real for v in pivots[2]]
    moments,heat,(field,)=lib("three-index-free-atoms")["broken_fields"](root,[number])
    ctx.prec=160
    assert ordered[-1]<1
    pos=sum(v for v in c if v>0)
    neg=sum(-v for v in c if v<0)
    assert all(v>0 or v<0 for v in c)
    floor=arb("-.00139")
    ceiling=arb(".000007")
    def certify(slope,stop,lower_bound):
        # 对 F=K+λW 有 F'=(4π)^{-1} sum_s (2α_s-λg_s x) exp(-z_s/x)。
        # 从区间的中点作均值展开。原点的指数按质量分别取连续单侧界。
        stack=[(arb(0),stop)]
        leaves=0
        while stack:
            lo,hi=stack.pop()
            mid=(lo+hi)/2
            x=lo.union(hi)
            deriv=arb(0)
            for g,a,z in field:
                if z.is_zero():
                    e=arb(1)
                elif lo.is_zero():
                    e=arb(0).union((-z/hi).exp())
                else:
                    e=(-z/lo).exp().union((-z/hi).exp())
                deriv+=(2*a-slope*g*x)*e
            k,w=moments(field,heat,mid)
            value=k+slope*w+arb(0,((hi-lo)*abs(deriv)/(8*arb.pi())).upper())
            if (value>floor if lower_bound else value<ceiling):
                leaves+=1
                continue
            assert hi-lo>arb("1e-9"),(lo,hi,value)
            stack.extend(((lo,mid),(mid,hi)))
        return leaves
    lam=arb("107.49")
    first=certify(lam,arb(".00205"),False),certify(lam,arb(1),True)
    second=certify(arb(0),arb(1),False),certify(arb(0),arb(".000515"),True)
    implied=pos*ceiling-neg*floor
    assert implied<1
    print("正权上界假设、负权上界假设的两套闭区间数；同一个有符号和上界",
          first,second,implied.upper().str(12),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
