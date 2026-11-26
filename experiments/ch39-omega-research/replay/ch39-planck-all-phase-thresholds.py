"""三重轮换谱带在两种普朗克标定下的全相位阈值排除。"""
from functools import cache
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


ROOT=Path(__file__).parent
phase=run_path(str(ROOT/"ch39-neutral-mean-phase-inertia.py"))
local=run_path(str(ROOT/"ch39-neutral-mean-local-mass.py"))
bands={
    1:"1/3.5/1 5/13/6 14.9891/14.9895/1 15.0015/15.00175/1 15.0027/15.0028/1 15.35/16/2 16/19/3 19/20/3 22.1/23.1/1",
    2:".04/3.5/1 5/13/6 15.1/15.14/1 15.25/15.35/1 15.7/18/6 18/20.3/4",
    3:"1/3.5/1 5/13/6 14.935/14.943/1 14.9956/14.9959/1 15.013/15.0152/1 15.45/16/2 16/19/3 19/20/3 22.1/23.1/1",
    4:".14/3.5/1 5/13/6 15.009/15.024/1 15.25/15.45/1 15.7/18/6 18/20.3/4",
    6:".33/3.5/1 5/13/6 14.89/14.91/1 15.35/15.5/1 15.7/18/6 18/21.9/4"}


def dimensionless(model):
    @cache
    def count(q,s):
        return phase["count_below"](q,arb(s))
    # 三条径向边的总算子扰动不超过 4×单边改变量。c 取实际总迹/18。
    offset=model["shell"]["soft"]["point"]["study"]["offset"]
    assert offset>0
    assert 4*arb(23)/25*arb(2)**-20/(14**2+offset) < arb("1e-7")
    pad=arb("1e-7")
    result=[]
    for q,table in bands.items():
        n=0
        for entry in table.split():
            a,b,k=entry.split("/")
            k=int(k)
            assert count(q,a)==n and count(q,b)==n+k
            n+=k
            left,right=arb(a)-pad,arb(b)+pad
            hs=[v*abs(v-15) for v in (left,right)]
            for turning in (arb("7.5"),arb(15)):
                if left < turning < right:
                    hs.append(turning*abs(turning-15))
            hull=hs[0]
            for h in hs[1:]:
                hull=hull.union(h)
            # 使用 |c-15| < pad，剩余宽度涵盖中心随零点改变。
            result.append((hull+arb(0,right*pad),q,k))
        assert n==19
    for v in model["shell"]["shell_spectra"]((0,))[0]:
        v=v+arb(0,2*pad)
        result.append((v*abs(v-15)+arb(0,v*pad),0,1))
    assert all(h>0 for h,_,_ in result)
    return result


def verify():
    model,_,_,initial,weights=local["inputs"]()
    ctx.prec=192
    hs=dimensionless(model)
    hmax=arb(0)
    for h,_,_ in hs:
        hmax=hmax.max(h)
    assert 156<hmax<188
    directions=((1,0,-1),(0,1,-1),(-1,1,0),(-1,0,1),(-5,-3,8),(-2,-3,5),(-3,5,-2))
    def dot(a,b):
        return sum(x*y for x,y in zip(a,b))
    # 对低能强耦合的一整个区间做检验。
    initial[2]=1/arb(".117").union(arb(".119"))
    mz,pi=model["mz"],arb.pi()
    logmin,logmax=[(arb(x)/mz).log() for x in ("1e13","1e17")]
    slices=256
    times=[logmin+(logmax-logmin)*i/slices for i in range(slices+1)]
    for normalize in (False,True):
        unit=arb("2.435e18")/(hmax if normalize else 1)
        rows=[((unit.lower()*h.lower()).union(unit.upper()*h.upper()),q,k) for h,q,k in hs]
        assert all(m>mz for m,_,_ in rows)
        worst=[None,None]
        def exposures(fn,m,scale):
            # 在各端点使用对数核对匹配动量递增、对平方质量递减。
            rlo=(m.lower()/(mz*scale.upper().exp()))**2
            rhi=(m.upper()/(mz*scale.lower().exp()))**2
            lo=fn(rhi,arb(1))
            hi=fn(rlo,arb(1))
            ref=mz**2/(8*m.lower()**2)
            return lo.union(hi)-arb(0).union(ref.upper())
        vec=lambda x,p:model["sheet"](x,p,True)
        fer=lambda x,p:model["sheet"](x,p,False)
        for a,b in zip(times,times[1:]):
            scale=a.union(b)
            sharp={}
            smooth={}
            for m,q,k in rows:
                sharp[q]=sharp.get(q,0)+k*(scale-(m/mz).log()).max(0)
                kind=vec if q==0 else fer
                smooth[q]=smooth.get(q,0)+k*exposures(kind,m,scale)
                if q==3:
                    smooth["H"]=smooth.get("H",0)+k*exposures(model["scalar"],m,scale)
            base=[2*pi*x-v*scale for x,v in zip(initial,model["beta"])]
            scores=[[],[]]
            for d in directions:
                projected=dot(d,base)
                boson=dot(d,model["higgs"])
                sv=projected-sum(dot(d,weights[q])*h for q,h in sharp.items())
                fv=projected-dot(d,weights[0])*smooth[0]-boson*smooth["H"]
                fv-=sum(6*dot(d,v)*smooth[q] for q,v in zip(model["charges"],model["weyl"]))
                for bag,v in zip(scores,(sv,fv)):
                    bag.append((v/(pi*sum(map(abs,d)))).lower())
            for i,bag in enumerate(scores):
                bound=max(bag)
                worst[i]=bound if worst[i] is None else min(worst[i],bound)
        print("重端归一",normalize,"尖阈值与有限减法的逆耦合极差下界",*[v.str(12) for v in worst])
        bounds=(arb("2.037"),arb("3.334")) if normalize else (arb("3.6039"),arb("3.6039"))
        assert all(v>bound for v,bound in zip(worst,bounds))


if __name__=="__main__":
    verify()
