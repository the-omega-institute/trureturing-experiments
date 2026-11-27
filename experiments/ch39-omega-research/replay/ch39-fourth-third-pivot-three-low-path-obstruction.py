"""第四壳原生第三枢轴至多三个低截止时，中性端点匹配强制中途色逆耦合为负。"""
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    n=2_306_888_100_085
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,n+4) if row[0]==6)
    c=[a.real for a in pivots[2]]
    assert sum(c).contains(1)
    _,ordered=lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=n+1)[2]
    gauge=lib("all-field-signed-crossing")
    neutral=lib("full-neutral-free-cutoff")
    data=gauge["joint_system"](root,"1.176397e14",number=n+4,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,n+4)
    ctx.prec=160
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    p0,p1=data["p0"],data["p1"]
    pm=(arb(5600)/arb("1.22089e19"))**2
    lo=5*p0/4
    assert ordered[0]>lo
    cut=arb("1e-6")
    pi=arb.pi()
    kappa,eta=arb(".172"),arb(".168")
    a1,a2,a3=data["initial"]
    starts=[5*a1/3+a2-8*a3/3,a2-a3,arb(0),a3]
    def project(r):
        p,m,d,k=r
        return [p,m,m-kappa*p,m-kappa*p-15*d,m-eta*p,d,k]
    rule=gauge["angular_rule"](11,32)
    def coeffs(b):
        qq,mix,col=b
        return [(qq-8*col/arb(3))/(4*pi),(mix-col)/(4*pi),arb(0),col/(4*pi)]
    coeff=[project(coeffs(b)) for _,_,b in species]
    series=native["weak_series"](56)
    def weak(x):
        return series(x/p0) if x<1e7*p0 else neutral["weak_column"](x,p0,p1)
    def sample(x):
        exposed=[]
        for kind,z,b in species:
            high=gauge["angular_exposure"](kind,z,x,p0,p1,(p1/p0).log(),rule)
            middle=gauge["angular_exposure"](kind,z,x,p0,pm,(pm/p0).log(),rule) if b[2] else arb(0)
            exposed.append((high,middle))
        return exposed
    nproj=7
    weak_coef=(0,0,0,-15,0,1,0)
    high_ref=sample(cut)
    dref=weak(cut)
    base=project(starts)
    for j in range(nproj):
        base[j]+=weak_coef[j]*dref-sum(b[j]*v[int(j==6)] for b,v in zip(coeff,high_ref,strict=True))
    # 高端各层共享参照截止 cut，偏移半径来自积分 dR/dlog(x)<=delta_p/(4x)。
    # D 的参考常数相消，保留两端的有限截止误差。
    residual=[sum(abs(b[j]) for b in coeff)*((pm if j==6 else p1)-p0)/(4*cut)
              +abs(weak_coef[j])*2*(arb("2.398e-21")+arb("7.468e-36")/cut) for j in range(nproj)]
    span=(cut/lo).log()
    @lru_cache(maxsize=None)
    def at(index):
        # index 是对数区间的二进有理数，公共端点只计算一次。
        x=lo*(span*index).exp()
        values=sample(x)
        raw=[-sum(b[j]*(v[int(j==6)]-h[int(j==6)]) for b,v,h in zip(coeff,values,high_ref,strict=True))
             for j in range(nproj)]
        return x,raw,values

    def moment_bound(p,z,left,right):
        # y exp(-y)<=min(1/e, y, 2exp(-y/2))，两个对称半区的指数角平均至多 24x/p。
        return min(arb(".368"),((z+p/4)/left).upper(),(2*(-z/(2*right)).exp()).upper(),(48*right/p).upper())

    @lru_cache(maxsize=None)
    def band(depth,index):
        den=2**depth
        left,gl,vl=at(index/den)
        right,gr,vr=at((index+1)/den)
        interval=left.union(right)
        # 跨越两个弱差表示的接缝时分别按各自的有效域取界。
        if right<1e7*p0:
            d=series(interval/p0)
        elif left>=1e7*p0:
            d=neutral["weak_column"](left,p0,p1)
        else:
            d=series(left.union(1e7*p0)/p0).union(neutral["weak_column"](1e7*p0,p0,p1))
        dw=d-dref
        result=[]
        sq=(span/den)**2/8
        for j in range(nproj):
            if j==5:
                result.append(dw)
                continue
            p=pm if j==6 else p1
            ix=int(j==6)
            # 两个端点的线性插值误差界，以及逐物种单调的直接界，取交集。
            error=sq*sum(abs(b[j])*(moment_bound(p0,z,left,right)+moment_bound(p,z,left,right))
                         for b,(_,z,_) in zip(coeff,species,strict=True))
            chord=gl[j].union(gr[j])+arb(0,error.upper())
            mono=-sum(b[j]*(v[ix].union(w[ix])-h[ix]) for b,v,w,h in zip(coeff,vl,vr,high_ref,strict=True))
            result.append(chord.intersection(mono)+weak_coef[j]*dw)
        return tuple(result)

    # 先用 64 个对数子区间生成粗树，之后只在尚可能为真的闭盒上继续二分。
    @lru_cache(maxsize=None)
    def bounds(depth,index):
        if depth>=6:
            return band(depth,index)
        return tuple(l.union(r) for l,r in zip(bounds(depth+1,2*index),bounds(depth+1,2*index+1),strict=True))

    counts=[0,0]
    max_depth=0
    # 少于三层时把缺位的坐标设为参照 cut。全体权的误差预算仍包括每个实际高层。
    total=sum(map(abs,c))
    for active in combinations(range(19),3):
        constants=[v+arb(0,(total*e).upper()) for v,e in zip(base,residual,strict=True)]
        stack=[((0,0),)*3]
        while stack:
            box=stack.pop()
            ranges=[bounds(*cell) for cell in box]
            counts[0]+=1
            vals=[v+sum(c[i]*r[j] for i,r in zip(active,ranges,strict=True)) for j,v in enumerate(constants)]
            if any(not v.contains(0) for v in vals[:-1]):
                continue
            if vals[-1]<-1:
                counts[1]+=1
                continue
            scores=[float(abs(c[i])*(r[0].rad()+20*r[2].rad()+150*r[5].rad())) for i,r in zip(active,ranges,strict=True)]
            choices=[k for k,(d,_) in enumerate(box) if d<19]
            assert choices,(tuple(i+1 for i in active),box,[v.str(8) for v in vals])
            k=max(choices,key=lambda k:scores[k])
            depth,index=box[k]
            max_depth=max(max_depth,depth+1)
            for side in (0,1):
                stack.append((*box[:k],(depth+1,2*index+side),*box[k+1:]))
    print("三行组合、总盒、色负排除盒、最大深度、角积分使用的不同截止端点数",
          969,*counts,max_depth,at.cache_info().currsize,flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
