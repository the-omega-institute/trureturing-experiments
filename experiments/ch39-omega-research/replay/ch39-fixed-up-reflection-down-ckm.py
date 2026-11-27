"""固定右侧上型反射，在全部上型中心质量三平面上检验任意下型原生谱权。"""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def fixed_moments(root):
    base = run_path(str(root / "ch39-rank-two-generator-up-quark.py"))
    training = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    transport = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    shells = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    orbits = run_path(str(root / "ch39-up-quark-arbitrary-site-phase-floor.py"))["ORBITS"]
    n = base["N"]
    _,labels,_ = training["verify_nodes"](root, ((n+1,".0294"),), ".1", (),
                                         ((6,7,10,13),(6,7,10,12)))
    ctx.prec = 384
    p = base["projector"](orbits)
    reflection = base["eye"](19)-2*p
    _,source = transport["spectrum"](shells,-4,n+1)
    _,target = transport["spectrum"](shells,-4,n+2)
    weights = base["read"](source,target,labels[n+1,-4],reflection)
    sites = [shells["soft"]["pole"],*shells["boundary_sites"]()]
    c = acb_mat([[target[i,j] for j in range(19)] for i in sites])
    raw = acb_mat([[c[i,j]*weights[j] for j in range(19)] for i in range(4)])*base["adj"](c)
    omega = acb(-arb(1)/2,arb(3).sqrt()/2)
    eigen = []
    for phase in (omega,omega.conjugate()):
        v = acb_mat([[0],[1],[phase],[phase**2]])*(1/arb(3).sqrt())
        z = (base["adj"](v)*raw*v)[0,0].real
        assert (raw*v-v*z).contains(acb_mat(4,1))
        eigen.append(z)
    h = arb(3).sqrt()
    trivial = acb_mat([[1,0],[0,1/h],[0,1/h],[0,1/h]])
    bb = base["adj"](trivial)*raw*trivial
    mid = (bb[0,0].real+bb[1,1].real)/2
    diff = (bb[0,0].real-bb[1,1].real)/2
    spread = (diff**2+abs(bb[0,1])**2).sqrt()
    eigen.extend((mid-spread,mid+spread))
    mass = list(map(arb,base["CENTERS"]))
    normal = [prod(a-m for m in mass)/prod(a-b for k,b in enumerate(eigen) if k!=i)
              for i,a in enumerate(eigen)]
    assert all(z>0 for z in normal) and sum(normal).contains(1)
    a,b = normal[:2]
    e,f = eigen[:2]
    mean = (b*e+a*f)/(a+b)
    variance = (1-a-b)*a*b*(e-f)**2/(a+b)**2
    assert mean > 0 and variance > 0
    print("四站本征值，依次为两字符、平凡块升序",[z.str(15) for z in eigen],flush=True)
    print("法向对应概率",[z.str(12) for z in normal],flush=True)
    print("共同第二弱轴的上型均值与标准差",mean.str(14),variance.sqrt().str(14),flush=True)
    return mean,variance.sqrt()


def verify(root):
    mean, deviation = fixed_moments(root)
    pencil = run_path(str(root / "ch39-triple-geometry-pencil.py"))
    lib,geom,v,mass,a,b,_ = pencil["setup"](root)
    for chart,(base,other) in enumerate(((a,b),(b,a))):
        stack=[(arb(-1),arb(1),(0,1,2),0)]
        count=depth=0
        least=arb(1)
        while stack:
            lo,hi,ids,level=stack.pop()
            count+=1
            depth=max(depth,level)
            lines=lib["lines"](base,other,lo,hi,ids,geom)
            keep=[]
            for i in ids:
                if i in lines:
                    s,_,variance=lib["moments"](v*lines[i],mass[0])
                    margin = abs(s-mean).lower()
                    if variance > 0:
                        margin = max(margin,abs(variance.sqrt()-deviation).lower())
                    if margin > arb("1e-5"):
                        least=min(least,margin)
                        continue
                keep.append(i)
            if keep:
                assert level<24,(chart,keep,lo,hi)
                mid=(lo+hi)/2
                stack.extend(((lo,mid,keep,level+1),(mid,hi,keep,level+1)))
        print("铅笔图、区间、最大二分深度、偏离量严格下界",chart,count,depth,least.str(15),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
