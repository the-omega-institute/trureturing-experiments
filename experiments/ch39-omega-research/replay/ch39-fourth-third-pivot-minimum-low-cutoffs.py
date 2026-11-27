"""第四壳实际第三枢轴在固定能标匹配三项中性读数所需的最少低截止数。"""
from functools import lru_cache
from math import factorial
from itertools import combinations
from pathlib import Path
from runpy import run_path
from flint import arb,ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number=2_306_888_100_089
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(row for row in native["shell_weights"](root,number) if row[0]==6)
    weights=[v.real for v in pivots[2]]
    assert sum(weights).contains(1)
    _,ordered=lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=number-3)[2]
    gauge=lib("all-field-signed-crossing")
    data=gauge["joint_system"](root,"1.176397e14",number=number,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec=192
    p0,p1=data["p0"],data["p1"]
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    neutral=lib("full-neutral-free-cutoff")
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    starts=data["initial"]
    startp=5*starts[0]/3+starts[1]-8*starts[2]/3
    startm=starts[1]-starts[2]
    rule=gauge["angular_rule"](16,43)
    weak_ref=arb("-.074256577194663365545")
    integrals=[]
    for kind,z,_ in species:
        val=sum(w[kind]*((z+p1*h)/(z+p0*h)).log() for h,w in rule[0])
        integrals.append(val+arb(0,(rule[1]*(p1/p0).log()).upper()))
    projections=((1,0),(0,1),(-arb(".172"),1))
    gauge_rule=gauge["angular_rule"](12,30)
    kernel=lib("fourth-free-first-probe-dual")
    constants=[weak_ref]
    gauge_at=[]
    for p,m in projections:
        coef=[-(p*(arb(q)-arb(8)*col/3)+m*(mix-col))/(4*arb.pi()) for _,_,(q,mix,col) in species]
        inf=sum(c*v for c,v in zip(coef,integrals,strict=True))
        constants.append(p*startp+m*startm+inf)
        at,size,_=kernel["gauge_combination"](species,gauge_rule,p0,p1,(p,m))
        gauge_at.append((at,size,inf))
    origin=arb("1.2768")*p0
    # 极低端的高动量 E1 项在所有 Gauss 节点上均可舍入到带正尾的零。
    # 对五个轻场在原 Gauss 节点积掉 E1 小自变量级数，避免大正负权抵消损失。
    exact_rule=gauge["angular_rule"](13,36)
    cutoff_switch=arb("1e-27")
    nodes,quaderr=exact_rule
    assert min(h*p1 for h,_ in nodes)>100*cutoff_switch
    assert all(z>100*cutoff_switch for _,z,_ in species[5:])
    exact=[]
    for p,m in projections:
        c=[-(p*(arb(q)-arb(8)*col/3)+m*(mix-col))/(4*arb.pi()) for _,_,(q,mix,col) in species]
        bound=sum(abs(q) for q in c)
        ae=[arb(0)]*23
        linear,constant=arb(0),arb(0)
        for zc,(kind,z,_) in zip(c[:5],species[:5],strict=True):
            for h,w in nodes:
                a=(z+p0*h)/p0
                b=zc*w[kind]
                linear+=b
                constant-=b*(arb.const_euler()+a.log())
                ap=arb(1)
                for n in range(1,23):
                    ap*=a
                    ae[n]+=b*ap
        ae=[ae[n]*(-1)**(n+1)/(n*factorial(n)) for n in range(1,23)]
        assert max((z+p0/4)/origin for _,z,_ in species[:5])<3
        err=bound*((p1/p0).log()*(quaderr+(-arb(100)).exp())+(-arb(100)).exp()+arb(3)**23/(23*factorial(23)))
        exact.append((constant,linear,ae,err))
    weak=native["weak_series"](56)
    floor=arb("1e-6")
    saturated=[neutral["weak_column"](floor,p0,p1)-weak_ref]+[arb(0,(size*(p1-p0)/(4*floor)).upper()) for _,size,_ in gauge_at]
    assert origin<ordered[0] and floor<ordered[-1]
    @lru_cache(None)
    def edge(i,power):
        return origin*(floor/origin)**(arb(i)/2**power)
    @lru_cache(None)
    def sample(i,power):
        if power and i%2==0:
            return sample(i//2,power-1)
        x=edge(i,power)
        if x<cutoff_switch:
            val=[]
            ratio=p0/x
            for (constant,linear,coefs,error),(_,_,inf) in zip(exact,gauge_at,strict=True):
                polys=arb(0)
                for c in reversed(coefs):
                    polys=(polys+c)*ratio
                val.append(constant+linear*(x/p0).log()+polys-inf+arb(0,error.upper()))
            return val
        return [at(x)-inf for at,_,inf in gauge_at]
    @lru_cache(None)
    def segment(i,depth):
        a,b=edge(i,depth),edge(i+1,depth)
        length=((floor/origin).log()/2**depth)
        weakband=neutral["weak_column"](a,p0,p1)-weak_ref
        if b<100000*p0 and length<arb(".1"):
            weakband=weak(a.union(b)/p0)-weak_ref
        left,right=sample(i,depth),sample(i+1,depth)
        ranges=[z.union(v)+arb(0,(size*length**2/(8*arb(1).exp())).upper())
                for z,v,(_,size,_) in zip(left,right,gauge_at,strict=True)]
        return [weakband,*ranges]
    # 每个自由层还可在高端；高端其余十七层共用零中心误差界。
    highmax=[abs(v).upper() for v in saturated]
    total_variation=sum(map(abs,weights))
    @lru_cache(None)
    def allseg(i,depth):
        if depth==-1:return saturated
        return segment(i,depth)
    def exclude_pair(i,j):
        c,d=weights[i],weights[j]
        bg=[v+arb(0,((total_variation-abs(c)-abs(d))*err).upper()) for v,err in zip(constants,highmax,strict=True)]
        # 一低/零低由 (低,高)、(高,高) 分支覆盖。
        work=[(0,k,0,l) for k in (-1,0) for l in (-1,0)]
        count=0
        while work:
            a,k,b,l=work.pop();count+=1
            first,second=allseg(a,k),allseg(b,l)
            if any(not (v+c*f+d*s).contains(0) for v,f,s in zip(bg,first,second,strict=True)):continue
            assert (k,l)!=(-1,-1) and max(k,l)<19,(i+1,j+1,a,k,b,l)
            # 平分贡献变化较大的一个对数区间，兼顾两层。
            left_weight=abs(float(c.mid()))**.35*(2.**-k if k>=0 else 0.)
            right_weight=abs(float(d.mid()))**.35*(2.**-l if l>=0 else 0.)
            if left_weight>=right_weight:
                work.extend((2*a+h,k+1,b,l) for h in (0,1))
            else:
                work.extend((a,k,2*b+h,l+1) for h in (0,1))
        return count
    pairs=list(combinations(range(19),2))
    checks=sum(exclude_pair(i,j) for i,j in pairs)
    assert len(pairs)==171
    print("已排除的候选对数、检查的闭矩形数",len(pairs),checks,flush=True)

    # 上界来自已认证的同一酉生成元，同时保持源谱、引力矩与完整味目标。
    witness=lib("fourth-third-pivot-joint")["verify"](root)
    for threshold in (floor,arb("9e-5")):
        below=[i+1 for i,x in enumerate(witness) if x<threshold]
        above=[i+1 for i,x in enumerate(witness) if x>threshold]
        assert below==[1,10,17] and len(above)==16
        print("阈值、见证中的低截止层、其余十六层的最小下界",
              threshold.str(6),below,min(witness[i-1].lower() for i in above).str(13),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
