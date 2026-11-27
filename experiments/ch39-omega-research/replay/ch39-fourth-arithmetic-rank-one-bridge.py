"""三条原生算术射线等权选取荷六两端谱线，并跨七维平凡空间生成秩一桥。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def select(vectors, indices, density):
    cols = {j:acb_mat([[vectors[i,j]] for i in range(19)]) for j in indices}
    scores = {j:(adj(v)*density*v)[0,0].real for j,v in cols.items()}
    best = max(indices,key=lambda j:scores[j].mid())
    gap = scores[best] - max(z.upper() for j,z in scores.items() if j!=best)
    assert gap>0
    return best, cols[best], scores[best], gap


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-"+name+".py")))

    ctx.prec = 384
    shell = lib("three-ray-chiral-boundary")
    sites = shell["geo"]["sites"]
    rotation = [sites.index((j,k,i)) for i,j,k in sites]
    four = (shell["soft"]["pole"],*shell["boundary_sites"]())
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)
    density = acb_mat(19,19)
    for coeff in gram[1:]:
        state = acb_mat([[coeff[four.index(i)] if i in four else 0] for i in range(19)])
        density += state*adj(state)/(adj(state)*state)[0,0].real/3
    assert density.trace().real.contains(1)

    n = 2_306_888_100_085
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    left,right = [spectrum(shell,6,i)[1] for i in (n+3,n+4)]
    group = lib("charge-six-projector-algebra-transport")["groups"]
    old,new = [group(v,rotation) for v in (left,right)]
    a,s,ws,ds = select(left,old[0],density)
    b,t,wt,dt = select(right,new[1]+new[2],density)
    assert (adj(s)*t)[0,0].contains(0)
    coeff = (adj(t)*density*s)[0,0]
    size = abs(coeff)
    assert size>0
    bridge = (t*adj(s))*coeff
    generator = bridge-adj(bridge)
    eye = lib("rank-two-generator-up-quark")["eye"](19)
    z = (eye+generator)*(eye-generator).inv()
    assert (adj(z)*z).contains(eye)
    trivial = acb_mat([[left[i,j] for j in old[0]] for i in range(19)])
    p = trivial*adj(trivial)
    # 四站支撑只含一个三点轨道，故 G 在各非平凡字符上的像至多一维。
    # Frobenius 范数在这个秩一跨块上等于算子范数，并对任何两端选线给出上界。
    r = acb_mat([[int(i==rotation[j]) for j in range(19)] for i in range(19)])
    p0 = (eye+r+r*r)/3
    omega = acb(-arb(1)/2,arb(3).sqrt()/2)
    plus = (eye+r*omega+r*r*omega**2)/3
    full = plus*density*p0
    cap = (adj(full)*full).trace().real.sqrt()
    assert cap>size and 6*cap.atan()<arb.pi()/2
    caps = [(2*k*cap.atan()).sin()**2 for k in (2,1,3)]
    print("任意两端本征线的桥强度上界与三个幂的概率上界",cap.str(16),*[v.str(15) for v in caps])
    angle = 2*size.atan()
    print("两端能序、目的字符、各自得分及领先余量",a+1,b+1,
          next(k for k in (1,2) if b in new[k]),*[v.str(15) for v in (ws,wt,ds,dt)])
    print("有向桥唯一非零奇异值与旋转角",size.str(15),angle.str(15))
    for k,op in ((-2,adj(z*z)),(1,z),(3,z**3)):
        off = (eye-p)*op*p
        prob = (k*angle).sin()**2
        assert (adj(off)*off).trace().real.overlaps(prob) and prob>0
        measured = abs((adj(t)*op*s)[0,0])**2
        assert measured.overlaps(prob)
        print("幂与最大转移概率",k,prob.str(16),flush=True)
    return z


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
