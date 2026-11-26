"""用整数四站法向和固定正有理谱权，正向复算六个夸克汤川与 CKM。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


INCREMENTS = {
    -4: {1:".016977439", 6:".011999003", 13:".024098414",
         16:".00001001904247", 19:"2.051869559", 17:".011055728"},
    2: {1:".030778455", 6:"116.2670305", 10:"1.494367523",
        16:".00068577302", 17:".010910881", 19:"2.00209590"},
}


def adj(m):
    return m.conjugate().transpose()


def orthonormal_frame():
    raw = (8,363,363,-697)
    length = arb(sum(i*i for i in raw)).sqrt()
    n = acb_mat([[i/length] for i in raw])
    e0 = acb_mat([[1],[0],[0],[0]])
    identity = acb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    reflection = identity-(n+e0)*adj(n+e0)*(1/(1+n[0,0]))
    frame = acb_mat([[reflection[i,j] for j in (1,2,3)] for i in range(4)])
    assert (acb_mat([list(raw)])*frame).contains(acb_mat(1,3))
    assert (adj(frame)*frame).contains(acb_mat([[int(i==j) for j in range(3)] for i in range(3)]))
    return frame


def verify(root):
    ray = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    flavor = run_path(str(root/"ch39-three-gram-source-flavor.py"))
    ctx.prec = 320
    height = acb.zeta_zero(2_306_888_100_085).imag
    frame = orthonormal_frame()
    sites = [ray["soft"]["pole"],*ray["boundary_sites"]()]
    assert sites == [14,8,7,6]
    outputs = []
    projectors = []
    for q, scale in ((-4, ".8616"), (2, ".01401")):
        shell = ray["laplacian"](q,height)
        eig, vectors = shell.eig(right=True,algorithm="rump")
        order = sorted(range(19),key=lambda j:eig[j].real.mid())
        assert all(v.imag.contains(0) for v in eig)
        assert all(eig[i].real<eig[j].real for i,j in zip(order,order[1:]))
        block = acb_mat([[arb("1e-9")*int(i==j) for j in range(3)] for i in range(3)])
        for layer, increment in INCREMENTS[q].items():
            index = order[layer-1]
            norm_sq = sum(abs(vectors[i,index])**2 for i in range(19))
            assert norm_sq>0 and arb(increment)>0
            v = adj(frame)*acb_mat([[vectors[i,index]] for i in sites])
            block += v*adj(v)*(arb(increment)/norm_sq)
        levels, p = flavor["projections"](block)
        outputs.extend(arb(scale)*value for value in levels)
        projectors.append(p)
    pu,pd=projectors
    probabilities=[[(a*b).trace().real for b in pd] for a in pu]
    s13=probabilities[0][2].sqrt()
    c13=(1-s13*s13).sqrt()
    s12,s23=probabilities[0][1].sqrt()/c13,probabilities[1][2].sqrt()/c13
    c12,c23=(1-s12**2).sqrt(),(1-s23**2).sqrt()
    j=(pu[0]*pd[0]*pu[1]*pd[1]).trace().imag
    sin_delta=j/(s12*s23*s13*c12*c23*c13**2)
    cos_delta=(probabilities[1][0]-(s12*c23)**2-(c12*s23*s13)**2)/(2*s12*s23*c12*c23*s13)
    assert sin_delta>0 and cos_delta>0 and (sin_delta**2+cos_delta**2).contains(1)
    delta=arb.atan2(sin_delta,cos_delta)
    angles = [s12.asin(),s23.asin(),s13.asin(),delta]
    published = (
        ("y_u",".00000615",".00000014"), ("y_c",".00311",".00005"), ("y_t",".8616",".0043"),
        ("y_d",".0000135",".0000002"), ("y_s",".000268",".000003"), ("y_b",".01401",".00009"),
        ("theta12",".22704",".00082"), ("theta23",".04275",".00042"),
        ("theta13",".00377",".00008"), ("delta","1.139",".023"),
    )
    deviations = []
    for (name,mean,radius),value in zip(published,outputs+angles,strict=True):
        center,halfwidth = arb(mean),arb(radius)
        assert center-halfwidth < value < center+halfwidth
        deviations.append(abs(value-center)/halfwidth)
        print(name,value.str(15))
    worst = max(d.upper() for d in deviations)
    assert worst<arb(".30")
    print("相对于所列边际半宽的最大偏差、有符号 J",worst.str(13),j.str(15))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
