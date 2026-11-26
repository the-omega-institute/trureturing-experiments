"""四个真实算术 Gram 源每次保留三个，对任意非负上下谱响应检验质量与顶底夹角。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path
from flint import acb, arb, arb_mat, ctx


def inner(v,w):
    return sum(a.conjugate()*b for a,b in zip(v,w,strict=True))


def normal(gram,omit):
    kept=[r for i,r in enumerate(gram) if i!=omit]
    vec=[(-1)**i*arb_mat([[v for j,v in enumerate(row) if j!=i] for row in kept]).det()
         for i in range(4)]
    length=sum(v*v for v in vec).sqrt()
    return [acb(v/length) for v in vec]


def project(vec,n):
    f=inner(n,vec)
    vals=[a-f*b for a,b in zip(vec,n,strict=True)]
    length=inner(vals,vals).real.sqrt()
    assert length>0
    return [v/length for v in vals]


def angle(v,w):
    s=sum(abs(v[i]*w[j]-v[j]*w[i])**2 for i,j in combinations(range(4),2))
    sin=arb(0).union(s.upper().sqrt()) if s.contains(0) else s.sqrt()
    return arb.atan2(sin,abs(inner(v,w)))


def tighten(rays,i,r,ratio):
    # 完全相同的射线可合并；极近但未证明相同的方向保留其最大角误差。
    delta=[angle(rays[i],v) for v in rays]
    near=[t for t in delta if t<arb("1e-5")]
    far=[t for t in delta if not t<arb("1e-5")]
    eps=max(t.upper() for t in near)
    if any(not t>r for t in far):
        return r
    for _ in range(8):
        gain=max((t.sin()/(t-r).sin()**2).upper() for t in far)
        predicted=(eps.sin()*(1+ratio)+ratio*gain).asin().upper()
        r=min(r,predicted)
    return r


def verify(root):
    ray=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    ctx.prec=320
    gram=native["gram_cholesky"](root)
    ckm=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ctx.prec=320
    def band(c,w):
        return arb(c)+arb(0,arb(w).upper())
    s23,s13=[v+arb(0,ckm["sine_radius"].upper()) for v in ckm["centers"][1:3]]
    book=((1-s23**2)*(1-s13**2)).sqrt().acos()
    table=(band(".04275",".00042").cos()*band(".00377",".00008").cos()).acos()
    for window in (book,table):
        assert arb(".0419483")<window<arb(".0434")
    bounds=((".04122939",".13307639"),(".04139032",".13476202"),
            (".04156149",".13672882"),(".04174777",".13897059"))
    sites=[ray["soft"]["pole"],*ray["boundary_sites"]()]
    height=acb.zeta_zero(2_306_888_100_085).imag
    raw=[]
    for q in (-4,2):
        vals,vec=ray["laplacian"](q,height).eig(right=True,algorithm="rump")
        seq=sorted(range(19),key=lambda i:vals[i].real.mid())
        assert all(vals[i].real<vals[j].real for i,j in zip(seq[:-1],seq[1:],strict=True))
        raw.append([[vec[k,i] for k in sites] for i in seq])
    ratio=[(arb(a)**2+arb(b)**2)/arb(c)**2 for a,b,c in
           ((".00000629",".00316",".8573"),(".0000137",".000271",".01392"))]
    radii=[(r/(1+r)).sqrt().asin().upper() for r in ratio]
    for omit in range(4):
        n=normal(gram,omit)
        up,down=[[project(z,n) for z in sector] for sector in raw]
        ru,rd=[[tighten(vv,i,r,s) for i in range(19)] for vv,r,s in
               zip((up,down),radii,ratio,strict=True)]
        small,large,pending=[],[],[]
        for i,u in enumerate(up):
            for j,d in enumerate(down):
                a=angle(u,d)
                lower=a.lower()-ru[i]-rd[j]
                upper=a.upper()+ru[i]+rd[j]
                if upper<arb(".0419483"):
                    small.append(upper)
                elif lower>arb(".0434"):
                    large.append(lower)
                else:
                    pending.append((i+1,j+1,float(lower),float(upper)))
        assert (len(small),len(large),len(pending))==(74,287,0)
        lo,hi=map(arb,bounds[omit])
        assert max(small)<lo and min(large)>hi
        print("省略源",omit,"过小、过大、余留",len(small),len(large),pending,
              "两支界",max(small).str(10),min(large).str(10),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
