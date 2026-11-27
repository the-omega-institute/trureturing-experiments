"""主对数给出的连续第三步在第四壳的十九站与四站质量下界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def multiply_columns(matrix,scales):
    return acb_mat([[matrix[i,j]*scales[j] for j in range(19)] for i in range(19)])


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    joint = lib("next-common-ckm-gravity")
    cx,rx,cy,ry,_ = joint["verify"](root)
    rank = lib("rank-two-generator-up-quark")
    n = rank["N"]
    emb = lib("common-embedding-central-right-ckm")
    _,labels,_ = emb["verify_nodes"](root,((n+1,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    coords = lib("independent-unitary-quark-right-holdout")["coordinates"]
    rows = lib("free-four-site-native-positive-spectrum")["native"](root,n+4)
    spec = lib("adjacent-pvm-transport")["spectrum"]
    shell = lib("three-ray-chiral-boundary")
    windows = lib("three-shell-outer-holdout")["WINDOWS"]
    local = lib("next-common-ckm")["local"]
    seed = lib("one-generator-gravity-up-strange-bottom")["cayley"]
    ctx.prec = 384
    eye = rank["eye"](19)
    def generator(center,radius,name,directions):
        box = [v+arb(0,radius.upper()) for v in center]
        base = seed(root,rank["eye"],name)
        return local(base,box,directions,eye,False)[0]
    x = generator(cx,rx,"ch39-full-ckm-gravity-cayley.tsv",lib("one-generator-gravity-full-ckm")["DIRECTIONS"])
    y = generator(cy,ry,"ch39-next-common-ckm-gravity-cayley.tsv",joint["DIRECTIONS"])
    ratio = y*adj(x)
    eig,raw = ratio.eig(right=True,algorithm="rump")
    gap = min(abs(1+v).lower() for v in eig)
    separation = min(abs(a-b).lower() for i,a in enumerate(eig) for b in eig[i+1:])
    assert gap>arb(".17335") and separation>0
    angles = [v.log().imag for v in eig]
    theta = max(abs(a).upper() for a in angles)
    assert theta<arb("2.969")
    lengths = [sum(abs(raw[i,j])**2 for i in range(19)).sqrt() for j in range(19)]
    beam = multiply_columns(raw,[1/z for z in lengths])
    assert (adj(beam)*beam).contains(eye)
    print("到主对数切口的最小单位弦距和最大绝对主辐角",gap.str(12),theta.str(12),flush=True)

    sectors = []
    for iq,(charge,projs) in enumerate(zip((-4,2),rows,strict=True)):
        source,dest = [spec(shell,charge,n+j)[1] for j in (1,4)]
        weights = labels[n+1,charge]
        spread = max(w.upper() for w in weights)-min(w.lower() for w in weights)
        if iq==0:
            left = adj(dest)*adj(y)*beam
            middle = adj(beam)*adj(y)*beam
            right = adj(beam)*adj(y*y)*adj(x*x)*source
        else:
            left = adj(dest)*beam
            middle = None
            right = adj(beam)*y*y*x*source
        cells = arb_mat([coords(p) for p in projs]).transpose()
        sectors.append((left,middle,right,weights,spread,cells))

    mesh = 1024
    for iq,(left,middle,right,weights,spread,cells) in enumerate(sectors):
        power = (2,1)[iq]
        lower = [arb(1)]*4
        upper = [arb(1)]*4
        raw_low = arb(1)
        accessible = {}
        for j in range(mesh+1):
            phases = [acb(0,(-1 if iq==0 else 1)*a*j/mesh).exp() for a in angles]
            at = multiply_columns(left,phases)
            if iq==0:
                at = multiply_columns(at*middle,phases)
            ov = at*right
            w = [sum(a*abs(ov[i,k])**2 for k,a in enumerate(weights)) for i in range(19)]
            raw_low = min(raw_low,min(v.lower() for v in w))
            trace,xx,yy,zz,p,m = (cells*arb_mat([[a] for a in w])).entries()
            r = (xx**2+yy**2+zz**2).sqrt()
            levels = [(trace-r)/2,(trace+r)/2,p,m]
            lows = sorted(v.lower() for v in levels)
            highs = sorted(v.upper() for v in levels)
            lower = [min(a,b) for a,b in zip(lower,lows,strict=True)]
            upper = [min(a,b) for a,b in zip(upper,highs,strict=True)]
            # 允许随质量选择三平面时，逐一检验该能级的压缩交错区间。
            for idx,(lo,hi) in enumerate(windows[iq]):
                if highs[idx]<arb(lo) and lows[idx+1]>arb(hi):
                    accessible.setdefault(idx,j)
        error = power**2*theta**2*spread/(4*mesh**2)
        nearest = power*theta*spread/(2*mesh)
        floor = [lower[0]-error]+[v-nearest for v in lower[1:]]
        global_floor = raw_low-error
        print("荷、整十九层最小下界、四站前两小本征值的连续下界与网格最小值上包",
              (-4,2)[iq],*[v.str(12) for v in (global_floor,floor[0],floor[1],upper[0],upper[1])],
              "第三质量可入整窗的节点",accessible.get(2),"/",mesh,flush=True)
        assert 2 in accessible
        if iq==0:
            assert global_floor>arb(".00118") and floor[0]>arb(".002187") and floor[1]>arb(".0225")
            assert upper[0]<arb(".002198") and upper[1]<arb(".02530")
        else:
            assert global_floor>arb(".000155") and floor[0]>arb(".000259") and floor[1]>arb(".00144")
            assert upper[0]<arb(".000261") and upper[1]<arb(".00188")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
