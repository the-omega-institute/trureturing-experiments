"""用荷六相邻本征线的极分解输运已校准的上一生成元。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(v):
    return v.conjugate().transpose()


def bridges(source, target, left, right):
    overlap = adj(target)*source
    weight = [[abs(overlap[j,i])**2 for j in range(19)] for i in range(19)]
    ordered,best = [],[]
    for old,new in zip(left,right,strict=True):
        ordered.extend(zip(old,new,strict=True))
        options = [(sum(weight[i][j] for i,j in zip(old,p,strict=True)),p) for p in permutations(new)]
        score,order = max(options,key=lambda x:x[0].mid())
        margin = score.lower()-max(v.upper() for v,p in options if p!=order)
        assert margin>0
        print("字符块维数与最大重叠配对优势",len(old),margin.lower().str(10),flush=True)
        best.extend(zip(old,order,strict=True))
    assert ordered == best
    for name,pairs in (("能序",ordered),("重叠",best)):
        d = acb_mat(19,19)
        for i,j in pairs:
            coeff = overlap[j,i]
            assert abs(coeff)>0
            d[j,i] = coeff/abs(coeff)
        w = target*d*adj(source)
        print(name,"19 对谱线的最小重叠模",min(abs(overlap[j,i]).lower() for i,j in pairs).str(10))
        yield name,w


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    joint = lib("next-common-ckm-gravity")
    cx,rx,cy,ry,_ = joint["verify"](root)
    rank = lib("rank-two-generator-up-quark")
    n,eye = rank["N"],rank["eye"](19)
    local = lib("next-common-ckm")["local"]
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"]
    def generator(name,c,r,directions):
        return local(seed(root,rank["eye"],name),[v+arb(0,r.upper()) for v in c],directions,eye,False)[0]
    x = generator("ch39-full-ckm-gravity-cayley.tsv",cx,rx,lib("one-generator-gravity-full-ckm")["DIRECTIONS"])
    y = generator("ch39-next-common-ckm-gravity-cayley.tsv",cy,ry,joint["DIRECTIONS"])
    _,labels,_ = lib("common-embedding-central-right-ckm")["verify_nodes"](
        root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    site = shell["geo"]["sites"]
    cycle = [site.index((b,c,a)) for a,b,c in site]
    groups = lib("charge-six-projector-algebra-transport")["groups"]
    spectra = {q:[spec(shell,q,v)[1] for v in (n+1,n+4)] for q in (-4,2,6)}
    src,dst = spectra[6]
    mid = spec(shell,6,n+3)[1]
    sz = [v.sqrt() for v in lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=n+1)[2][1]]
    moments,heat,fields = lib("three-index-free-atoms")["broken_fields"](root,[n+4])
    native = lib("native-rh-ray-pvm")
    _,_,pivots,_ = next(row for row in native["shell_weights"](root,n+4) if row[0]==6)
    c = [v.real for v in pivots[0]]
    ctx.prec = 384
    p0,p1 = [(arb(v)/arb("1.22089e19"))**2 for v in ("91.1876","1.176397e14")]
    weak = lib("full-neutral-free-cutoff")["weak_column"]
    windows = lib("three-shell-outer-holdout")["WINDOWS"]
    sites = (shell["soft"]["pole"],*shell["boundary_sites"]())
    assert sites==(14,8,7,6)
    for name,w in bridges(mid,dst,groups(mid,cycle),groups(dst,cycle)):
        z = w*y*adj(w)
        for q,op,limits in zip((-4,2),(adj(z*z)*adj(y*y)*adj(x*x),z*y*x),windows,strict=True):
            s,t = spectra[q]
            amps = rank["read"](s,t,labels[n+1,q],op)
            four = acb_mat([[t[i,j] for j in range(19)] for i in sites])
            raw = acb_mat([[four[i,j]*amps[j] for j in range(19)] for i in range(4)])*adj(four)
            ev = sorted((a.real for a in raw.eig(algorithm="rump")),key=lambda a:a.mid())
            assert all(a<b for a,b in zip(ev,ev[1:]))
            outside = [ev[i]>arb(hi) or ev[i+1]<arb(lo) for i,(lo,hi) in enumerate(limits)]
            print(name,q,"四站谱与交错排除",[a.str(12) for a in ev],outside,flush=True)
            assert outside == [True, True, False]
        f = rank["read"](src,dst,sz,z**3*y**3*x**3)
        assert min(f)**2>p0
        thermal = [moments(fields[0],heat,v*v) for v in f]
        columns = (*zip(*thermal),[weak(v*v,p0,p1) for v in f])
        vals = [sum(w*t for w,t in zip(c,col,strict=True)) for col in columns]
        print(name,"荷六 K、W、D",*[v.str(15) for v in vals],flush=True)
        assert all(not v.overlaps(goal) for v,goal in zip(vals,map(arb,("1","1e-122","0")),strict=True))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
