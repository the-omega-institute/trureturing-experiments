"""以选定中间能标的必要正性条件检验六个低截止的所有支撑。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx
import numpy as np


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
    probes=[(1334,"y"),(1334,"col"),(2371,"y"),(750000,"y"),(750000,"col"),
            (17800000000,"y"),(17800000000,"col"),(5600,"col"),(1800000,"col"),
            (56000000000,"y"),(10000000,"y"),(4200,"y"),(5600000,"col"),
            (100000000000,"y")]
    energies=list(dict.fromkeys(energy for energy,_ in probes))
    powers=[p1,*[(arb(v)/mp)**2 for v in energies]]
    lo=(arb("103.04188")/mp)**2
    cut=arb("1e-6")
    assert ordered[0]>lo
    span=(cut/lo).log()
    a1,a2,a3=data["initial"]
    def axes(rp,rm,w,col,y):
        return [rp,rm,rm-arb(".172")*rp,rm-arb(".172")*rp-15*w,rm-arb(".168")*rp,w,
                *[y if channel=="y" else col for _,channel in probes]]
    count=6+len(probes)
    beta=[axes((q-8*cc/arb(3))/(4*arb.pi()),(m-cc)/(4*arb.pi()),arb(0),
               cc/(4*arb.pi()),(q-m)/(4*arb.pi())) for _,_,(q,m,cc) in species]
    const=axes(5*a1/3+a2-8*a3/3,a2-a3,arb(0),a3,5*a1/3)
    which=[0]*6+[energies.index(energy)+1 for energy,_ in probes]
    weakcoef=[0,0,0,-15,0,1]+[0]*len(probes)
    kappa=arb("7.468e-36")
    tv=sum(map(abs,c))
    negative_mass=sum(-w for w in c if w<0)
    lib("neutral-positive-profile")["verify"]()
    # 负固有时密度给 D_inf <= D(x) <= 0，再按权号估计 D 的上端。
    # 对较高探测能标另用 D(M0)=0 将弱差改写为两个高能端点之差。
    weak_allowance=[]
    for j in (6+i for i,(_,channel) in enumerate(probes) if channel=="y"):
        pm=powers[which[j]]
        ent=lib("full-neutral-running")["weak_bound"]({"p0":pm},[pm],[arb(0)],
                                    lib("neutral-axial-top"),lib("neutral-bosonic"))
        a=(arb(".074256577194663365546")+ent)*negative_mass
        b=tv*(ent+2*(kappa/pm).sqrt())
        allowance=min(a.upper(),b.upper())
        const[j]+=allowance
        weak_allowance.append(allowance)
    # 浮点坐标和角权从 Arb 区间转出并逐项核对，线性和误差按最不利的顺序求和。
    nodes,qerr=gauge["angular_rule"](13,34)
    ev,table_error=lib("fourth-third-pivot-five-low-path-obstruction")["e1_interpolator"]()
    denom=[[] for _ in range(len(powers)+1)]
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
    drifts=[]
    for j in range(count):
        base.append(const[j]+ref[j]+weakcoef[j]*dref)
        drifts.append(massnorm[j]*(powers[which[j]]-p0)/(4*cut)+
                      abs(weakcoef[j])*2*(arb("2.398e-21")+kappa/cut))
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
                dw=wseries(x/p0).union(wseries(y/p0))-dref
            elif x>=limit:
                dw=neutral["weak_column"](x,p0,p1)-dref
            else:
                dw=wseries(x/p0).union(wseries(limit/p0)).union(neutral["weak_column"](limit,p0,p1))-dref
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
    high_error=np.array([[bracket(abs(w)*d)[1] for d in drifts] for w in c],dtype=np.int64)
    assert int(np.max(np.abs(high_error)))*19<2**50
    bounds=np.zeros((19,2*fine,count,2),dtype=np.int64)
    for i,w in enumerate(c):
        bounds[i,fine:]=[[bracket(w*v) for v in band] for band in graph]
    assert int(np.max(np.abs(bounds)))*7+int(np.max(np.abs(high_error)))*19+int(np.max(np.abs(constants)))<2**60
    for n in range(fine-1,0,-1):
        bounds[:,n,:,0]=np.minimum(bounds[:,2*n,:,0],bounds[:,2*n+1,:,0])
        bounds[:,n,:,1]=np.maximum(bounds[:,2*n,:,1],bounds[:,2*n+1,:,1])
    # 在临时目录编译深度优先遍历，内存不随同一层的未决盒数增长。
    from tempfile import TemporaryDirectory
    from subprocess import run
    import shutil
    compiler=shutil.which("c++")
    assert compiler,"六层整数覆盖需要 C++ 编译器"
    assert int(np.max(np.abs(bounds)))*400 < 2**63
    print("中性探测能标的弱差上界",[v.str(12) for v in weak_allowance],flush=True)
    with TemporaryDirectory() as folder:
        binary=Path(folder)/"cover"
        run([compiler,"-std=c++17","-O3",str(root/"ch39-low-cutoff-integer-cover.cpp"),"-o",str(binary)],check=True)
        with (Path(folder)/"table").open("w+b") as fp:
            np.array([6,fine,count],dtype=np.int64).tofile(fp)
            constants.tofile(fp)
            high_error.tofile(fp)
            bounds.tofile(fp)
            fp.seek(0)
            run([str(binary)],stdin=fp,check=True)



if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
