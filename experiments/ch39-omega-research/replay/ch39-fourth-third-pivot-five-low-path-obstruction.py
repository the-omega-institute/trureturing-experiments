"""以三个中间能标的必要正性条件检验五个低截止的所有支撑。"""
from itertools import combinations
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
import numpy as np


def e1_interpolator():
    # 二进区间内在 y 上作线性插值，E1''(y)=(1+y)e^-y/y² <=1/y²。
    div=512
    samples=[]
    for power in range(-205,6):
        edge=arb(2)**power
        for j in range(div):
            y=edge*(div+j)/div
            v=y.expint(1)
            f=float(v.mid())
            assert abs(v-f)<arb("1e-12")
            samples.append(f)
    samples.append(0.)
    table=np.array(samples)
    assert (-arb(64)).exp()/64<arb("1e-12")
    def ev(y):
        assert np.all(np.isfinite(y)) and np.min(y)>2.**-205
        m,k=np.frexp(np.minimum(y,64.))
        pos=(2*m-1)*div
        index=pos.astype(np.int64)
        frac=pos-index
        dest=(k.astype(np.int64)-1+205)*div+index
        return table[dest]+frac*(table[np.minimum(dest+1,len(table)-1)]-table[dest])
    return ev,arb(1)/(8*div**2)+arb("3e-11")


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    number=2_306_888_100_089
    native=lib("native-rh-ray-pvm")
    _,_,pivots,_=next(r for r in native["shell_weights"](root,number) if r[0]==6)
    c=[z.real for z in pivots[2]]
    assert sum(c).contains(1)
    _,ordered=lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=number-3)[2]
    gauge=lib("all-field-signed-crossing")
    neutral=lib("full-neutral-free-cutoff")
    data=gauge["joint_system"](root,"1.176397e14",number=number,order=4,last=8)
    signs,masses=lib("rotation-portal-top-energy")["rotation_signs"](root,number)
    ctx.prec=160
    mp=arb("1.22089e19")
    rho=(arb("246.22")/mp)**2/2
    species,_=neutral["response_rows"](data,masses,signs,rho,arb(".42995401872"))
    p0,p1=data["p0"],data["p1"]
    powers=[p1,(arb(5600)/mp)**2,(arb("1e7")/mp)**2,(arb("1e10")/mp)**2]
    lo=5*p0/4
    cut=arb("1e-6")
    assert ordered[0]>lo
    span=(cut/lo).log()
    a1,a2,a3=data["initial"]
    def axes(rp,rm,w,c1,hy,c2):
        return [rp,rm,rm-arb(".172")*rp,rm-arb(".172")*rp-15*w,rm-arb(".168")*rp,w,c1,hy,c2]
    count=9
    beta=[axes((q-8*cc/arb(3))/(4*arb.pi()),(m-cc)/(4*arb.pi()),arb(0),
               cc/(4*arb.pi()),(q-m)/(4*arb.pi()),cc/(4*arb.pi())) for _,_,(q,m,cc) in species]
    const=axes(5*a1/3+a2-8*a3/3,a2-a3,arb(0),a3,5*a1/3,a3)
    which=[0,0,0,0,0,0,1,2,3]
    weakcoef=(0,0,0,-15,0,1,0,0,0)
    # D(M0)=0。用 4x/p 和 H(p)+k/x 的较小者界定 D(P)-D(M0) 的单层读数。
    pm=powers[2];kappa=arb("7.468e-36")
    ent=lib("full-neutral-running")["weak_bound"]({"p0":pm},[pm],[arb(0)],
                                lib("neutral-axial-top"),lib("neutral-bosonic"))
    weak_mid=sum(map(abs,c))*(ent+2*(kappa/pm).sqrt())
    assert weak_mid<arb(".1")
    const[7]+=weak_mid.upper()

    # 浮点坐标和角权从 Arb 区间转出并逐项核对，线性和误差按最不利的顺序求和。
    nodes,qerr=gauge["angular_rule"](13,34)
    ev,table_error=e1_interpolator()
    denom=[[] for _ in range(5)]
    rows=[[] for _ in range(count)]
    coeff_error=[arb(0)]*count
    massnorm=[sum(abs(b[j]) for b in beta) for j in range(count)]
    for (kind,z,_),b in zip(species,beta,strict=True):
        for h,omega in nodes:
            for r,p in enumerate([p0,*powers]):
                val=z+p*h
                approx=float(val.mid())
                assert abs(approx/val-1)<arb("1e-13")
                denom[r].append(approx)
            for j in range(count):
                factor=-b[j]*omega[kind]
                vf=float(factor.mid())
                coeff_error[j]+=abs(factor-vf)
                rows[j].append(vf)
    den=np.array(denom);coeff=np.array(rows)
    norm_angular=[sum(abs(b[j]) for b in beta) for j in range(count)]
    eps=arb(2)**-52
    ne=den.shape[1]+6
    assert ne*eps<arb("1e-6")
    # 正角权的 Gauss 权和逼近一，直接对浮点系数绝对和作认证，避免依赖精确归一。
    norm_float=[]
    for j in range(count):
        s=float(np.abs(coeff[j]).sum())
        norm_float.append(arb(s)/(1-ne*eps))
    errs=[]
    for j in range(count):
        log=(powers[which[j]]/p0).log()
        round_sum=ne*eps/(1-ne*eps)*norm_float[j]*(log+1)
        errs.append((massnorm[j]*qerr*log+coeff_error[j]*(log+1)
                    +2*norm_float[j]*(table_error+arb("1e-12"))+round_sum).upper())
    def gauge_values(x):
        xf=float(x.mid())
        assert abs(xf/x-1)<arb("1e-13")
        E=ev(den/xf)
        deltas=E[0]-E[1:]
        return [arb(float(np.sum(coeff[j]*deltas[which[j]])))+arb(0,errs[j]) for j in range(count)]

    coarse,fine=512,4096
    assert fine%coarse==0 and fine&(fine-1)==0
    ref=gauge_values(cut)
    dref=neutral["weak_column"](cut,p0,p1)
    tv=sum(map(abs,c))
    wseries=native["weak_series"](56)
    base=[]
    for j in range(count):
        v=const[j]+ref[j]+weakcoef[j]*dref
        drift=tv*(massnorm[j]*(powers[which[j]]-p0)/(4*cut)+
                  abs(weakcoef[j])*2*(arb("2.398e-21")+kappa/cut))
        base.append(v+arb(0,drift.upper()))
    def end(j):
        return lo*(span*j/coarse).exp()
    def second(p,z,l,r):
        return min(arb(".368"),((z+p/4)/l).upper(),(2*(-z/(2*r)).exp()).upper(),(48*r/p).upper())
    graph=[]
    left=lo
    gl=gauge_values(left)
    for t in range(coarse):
        right=end(t+1);gr=gauge_values(right)
        bend=[]
        for j in range(count):
            p=powers[which[j]]
            e=sum(abs(b[j])*(second(p0,z,left,right)+second(p,z,left,right)) for b,(_,z,_) in zip(beta,species,strict=True))
            bend.append((e*(span/coarse)**2/8).upper())
        for j in range(fine//coarse):
            a=arb(j)/(fine//coarse);b=arb(j+1)/(fine//coarse)
            x=lo*(span*(t+a)/coarse).exp();y=lo*(span*(t+b)/coarse).exp()
            limit=10**7*p0
            if y<limit:
                dw=wseries(x.union(y)/p0)-dref
            elif x>=limit:
                dw=neutral["weak_column"](x,p0,p1)-dref
            else:
                dw=wseries(x.union(limit)/p0).union(neutral["weak_column"](limit,p0,p1))-dref
            graph.append([((1-a)*u+a*v).union((1-b)*u+b*v)-r+arb(0,e)+wc*dw
                          for u,v,r,e,wc in zip(gl,gr,ref,bend,weakcoef,strict=True)])
        left,gl=right,gr

    # 每个端点按 2^-32 向外取整。树中只做 int64 的最小、最大和加法。
    scale=2**32
    def bracket(v):
        lo=(v.lower()*scale).floor().unique_fmpz()
        hi=(v.upper()*scale).ceil().unique_fmpz()
        assert lo is not None and hi is not None
        return int(lo),int(hi)
    constants=np.array([bracket(v) for v in base],dtype=np.int64)
    bounds=np.zeros((19,2*fine,count,2),dtype=np.int64)
    for i,w in enumerate(c):
        bounds[i,fine:]=[[bracket(w*v) for v in band] for band in graph]
    assert int(np.max(np.abs(bounds)))*6+int(np.max(np.abs(constants)))<2**60
    for n in range(fine-1,0,-1):
        bounds[:,n,:,0]=np.minimum(bounds[:,2*n,:,0],bounds[:,2*n+1,:,0])
        bounds[:,n,:,1]=np.maximum(bounds[:,2*n,:,1],bounds[:,2*n+1,:,1])
    # 支撑与它的笛卡尔积盒按批处理。浮点分数只选择细分坐标。
    scores=np.zeros((19,2*fine),dtype=np.float64)
    for j,m in [(0,1),(2,20),(5,150)]:
        scores+=m*(bounds[:,:,j,1]-bounds[:,:,j,0])
    scores[:,fine:]=-1
    seen=excluded=0
    all_sets=np.array(list(combinations(range(19),5)),dtype=np.int64)
    for t in range(0,len(all_sets),128):
        active=all_sets[t:t+128]
        nodes=np.ones_like(active)
        while len(nodes):
            seen+=len(nodes)
            val=constants+bounds[active,nodes].sum(axis=1)
            keep=np.all(val[:,:6,0]<=0,axis=1)&np.all(val[:,:6,1]>=0,axis=1)
            pos=np.all(val[:,6:,1]>=0,axis=1)
            excluded+=int(np.count_nonzero(keep&~pos))
            active,nodes=active[keep&pos],nodes[keep&pos]
            if not len(nodes):
                break
            widths=scores[active,nodes]
            i=np.arange(len(nodes))
            j=widths.argmax(axis=1)
            bad=widths[i,j]<0
            assert not np.any(bad),(active[bad][0]+1,nodes[bad][0]-fine)
            nodes[i,j]*=2
            r=nodes.copy();r[i,j]+=1
            nodes=np.concatenate([nodes,r]);active=np.concatenate([active,active])
    print("五层支撑、全部闭盒、途中正性严格排除数、10^7 GeV弱差绝对界",
          comb(19,5),seen,excluded,weak_mid.upper().str(12),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
