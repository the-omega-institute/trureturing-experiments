"""用同一第三步精确保留第三枢轴荷六五式并实现三个下型质量中心。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


PARTITION = (1,1,1,1,1,0,1,1,1,1,2,1,2,1,2,1,2,0,1)


def adj(a):
    return a.conjugate().transpose()


def canon(frame):
    # 在每一谱线上选取第一个绝对值超过 0.02 的分量，将其相位取为正实。
    result=acb_mat(frame)
    for j in range(frame.ncols()):
        i=next(i for i in range(frame.nrows()) if abs(frame[i,j])>arb(".02"))
        assert all(abs(frame[k,j])<arb(".02") for k in range(i))
        phase=frame[i,j].conjugate()/abs(frame[i,j])
        for k in range(frame.nrows()):
            result[k,j]=frame[k,j]*phase
    return result


def leakage_bound(seed,shell,root_distance):
    sites=shell["geo"]["sites"]
    cycles=[sorted({i,sites.index((b,c,a)),sites.index((c,a,b))}) for i,(a,b,c) in enumerate(sites)]
    blocks=[orbit for i,orbit in enumerate(cycles) if i==orbit[0]]
    fixed=acb_mat([[1/arb(len(orbit)).sqrt() if i in orbit else 0 for orbit in blocks] for i in range(19)])
    assert len(blocks)==7
    out=seed*fixed
    off=out-fixed*(adj(fixed)*out)
    spectrum=(adj(off)*off).eig(algorithm="rump")
    upper=max(v.real.upper() for v in spectrum)
    return (upper.sqrt()+root_distance)**2


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    cyclic=lib("fourth-third-pivot-cyclic-five")
    amps,colour,x,y,incoming,readout,(low,h,k)=cyclic["verify"](root,return_inputs=True)
    ctx.prec=256
    inc=canon(incoming)
    target=canon(adj(readout))
    prob=[(a-low)/(h+k-2*low) for a in amps]
    assert all(q>0 for q in prob) and sum(prob).contains(1)
    p,q,r=[sum(v for v,i in zip(prob,PARTITION,strict=True) if i==j) for j in range(3)]
    cs=(r*r-p*p-q*q)/(2*p*q)
    assert -1<cs<1
    phase=acb(cs,(1-cs**2).sqrt())
    phases=[acb(1),phase,-(p+q*phase)/r]
    goal=target*acb_mat([[v.sqrt(),v.sqrt()*phases[i]] for i,v in zip(PARTITION,prob,strict=True)])
    assert (adj(goal)*goal).contains(acb_mat([[1,0],[0,1]]))
    rank=lib("rank-two-generator-up-quark")
    seed=lib("one-generator-gravity-up-strange-bottom")["cayley"](
        root,rank["eye"],"ch39-fourth-third-pivot-down-cayley.tsv")
    error=sum(abs(v).upper()**2 for v in (seed**3*inc-goal).entries()).sqrt().upper()
    # 用近恒等酉延拓将这对正交标架精确对齐，扰动 VZ0^3 的 Frobenius 范数。
    correction=2*arb(2).sqrt()*error
    levels=seed.eig(algorithm="rump")
    sep=min(abs(a*a+a*b+b*b).lower() for a in levels for b in levels)
    assert sep>arb(".0144")
    radius=arb(".0004")
    norm=(correction+3*radius**2+radius**3)/sep
    lip=(6*radius+3*radius**2)/sep
    print("两标架误差、立方导数最小奇异值、根映射半径和压缩比",
          *[v.str(15) for v in (error,sep,norm,lip)],flush=True)
    # 根 与 (根*)^-1 都在唯一性球中，因而该精确根是酉的。
    assert lip<1 and norm/(1-norm)<radius and 0<norm<1
    number=2_306_888_100_085
    _,labels,_=lib("common-embedding-central-right-ckm")["verify_nodes"](
        root,((number+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec=256
    shell=lib("three-ray-chiral-boundary")
    leakage=leakage_bound(seed,shell,norm)
    assert leakage<arb(".474")
    print("此联合解的最大跨空间转移概率上界",leakage.upper().str(15),flush=True)
    spec=lib("adjacent-pvm-transport")["spectrum"]
    source,dest=[spec(shell,2,number+j)[1] for j in (1,4)]
    weights=labels[number+1,2]
    assert all(w>=0 for w in weights)
    output=rank["read"](source,dest,weights,seed*y*x)
    four=acb_mat([[dest[i,j] for j in range(19)] for i in (14,8,7,6)])
    raw=lib("common-unitary-quark-right-trace")["weighted"](four,output)
    ev=sorted((v.real for v in raw.eig(algorithm="rump")),key=lambda v:v.mid())
    delta=norm*max(v.upper() for v in weights).sqrt()
    # 正的四站二次型由振幅线性算子的平方范数给出。其各奇异值都至多变化 delta。
    bands=[(v.sqrt()+arb(0,delta.upper()))**2 for v in ev]
    masses=list(map(arb,(".0000135",".000268",".01401")))
    assert all(bands[i]<m<bands[i+1] for i,m in enumerate(masses))
    print("四站谱上下界",[(v.lower().str(13),v.upper().str(13)) for v in bands],flush=True)
    print("荷六五式的色逆耦合区间",colour.lower().str(13),colour.upper().str(13),flush=True)
    return bands


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
