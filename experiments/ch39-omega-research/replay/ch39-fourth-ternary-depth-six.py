"""在 N+4 对只含 X,Y 叶子的无原点三元式穷举深度至六的全部指数。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def adj(a):
    return a.conjugate().transpose()


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    ks={-1,0}
    for d in range(1,7):
        differences={i-j for i in ks for j in ks}
        ks={m+k for m in differences for k in ks}
        assert ks==set(range(-(3**d+1)//2,(3**d-1)//2+1))
    assert len(ks)==730
    joint=lib("next-common-ckm-gravity")
    cx,rx,cy,ry,_=joint["verify"](root)
    rank=lib("rank-two-generator-up-quark")
    n=rank["N"]
    _,labels,_=lib("common-embedding-central-right-ckm")["verify_nodes"](
        root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    rows=lib("free-four-site-native-positive-spectrum")["native"](root,n+4)
    coords=lib("independent-unitary-quark-right-holdout")["coordinates"]
    ctx.prec=384
    eye=rank["eye"](19)
    local=lib("next-common-ckm")["local"]
    seed=lib("one-generator-gravity-up-strange-bottom")["cayley"]
    def generator(center,radius,dirs,file):
        z,_=local(seed(root,rank["eye"],file),center,dirs,eye,False)
        widths=[abs((v+arb(0,radius.upper()))-v).upper() for v in center]
        error=2*sum((1 if i==j else 2)*w*w for w,(i,j,_) in zip(widths,dirs,strict=True)).sqrt()
        return z,error.upper()
    x,dx=generator(cx,rx,lib("one-generator-gravity-full-ckm")["DIRECTIONS"],"ch39-full-ckm-gravity-cayley.tsv")
    y,dy=generator(cy,ry,joint["DIRECTIONS"],"ch39-next-common-ckm-gravity-cayley.tsv")
    ratio=y*adj(x)
    windows=lib("three-shell-outer-holdout")["WINDOWS"]
    spectrum=lib("adjacent-pvm-transport")["spectrum"]
    shell=lib("three-ray-chiral-boundary")
    sectors=[]
    for q,power,prior,projectors,win in zip((-4,2),(2,1),(adj(y*y)*adj(x*x),y*x),rows,windows,strict=True):
        old,dest=[spectrum(shell,q,index)[1] for index in (n+1,n+4)]
        w=labels[n+1,q]
        spread=max(v.upper() for v in w)-min(v.lower() for v in w)
        cells=arb_mat([coords(v) for v in projectors]).transpose()
        sectors.append((power,prior*old,adj(dest),w,spread,cells,win))
    names=("u","c","t","d","s","b")
    accepted={q:[] for q in names}
    floors={q:[] for q in names}
    for k in sorted(ks):
        z=(ratio**k if k>=0 else adj(ratio)**(-k))*y
        dz=abs(k)*(dx+dy)+dy
        for j,(power,src,read,w,spread,cells,win) in enumerate(sectors):
            op=adj(z*z) if power==2 else z
            ov=read*op*src
            out=arb_mat([[sum(a*abs(ov[i,l])**2 for l,a in enumerate(w))] for i in range(19)])
            tr,x1,x2,x3,p,m=(cells*out).entries()
            gap=(x1*x1+x2*x2+x3*x3).sqrt()
            levels=[(tr-gap)/2,(tr+gap)/2,p,m]
            wiggle=arb(0,(power*spread*(dx+dy+dz)).upper())
            lower=sorted((v+wiggle).lower() for v in levels)
            upper=sorted((v+wiggle).upper() for v in levels)
            for i,(lo,hi) in enumerate(win):
                name=names[3*j+i]
                floors[name].append((lower[i],k))
                if lower[i]>arb(hi) or upper[i+1]<arb(lo):
                    continue
                assert upper[i]<arb(hi) and lower[i+1]>arb(lo),(k,name)
                accepted[name].append(k)
    for q in names:
        inside=accepted[q]
        outside=sorted(ks-set(inside))
        value,index=min(floors[q])
        print("质量、排除数、可达数、最小的下端点及指数",q,len(outside),len(inside),value.str(13),index,flush=True)
        print("较短名单",("可达",inside) if len(inside)<=len(outside) else ("排除",outside),flush=True)
    return accepted


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
