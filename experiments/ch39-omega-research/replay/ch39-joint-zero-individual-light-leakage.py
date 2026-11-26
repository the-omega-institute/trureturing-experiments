"""固定十九站的对应关系，在每站任意调相时逐态认证轻空间的重泄漏。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


GROUPS = ((0,3,4), (1,2,5), (6,7,8), (9,15,16), (10,13,17), (11,12,18), (14,))
# 正对角 Cauchy–Schwarz 证书。按七个坐标轮换轨道用相同的正数。
CERTIFICATES = {
    2: (("20.475853","19.296504","13.663853",".422936",".423657","1.155762","16089.454604"),
        ".1687830",".171839"),
    -4: (("18.279364","17.378532","7.784257",".425672",".428529","1.165558","16221.902739"),
         ".1554432",".155712"),
    6: (("18.339628","17.660049","5.545251",".441459",".445593","1.163399","15266.831498"),
        ".148047023",".148047025"),
}


def compression_upper(basis, weights):
    scaled = acb_mat([[weights[i]*basis[i,j] for j in range(3)] for i in range(19)])
    compressed = basis.conjugate().transpose()*scaled
    # Gershgorin 行和对任意复叠加都给出最大本征值上界。
    return max(compressed[i,i].real.upper()
               +sum(abs(compressed[i,j]).upper() for j in range(3) if j!=i)
               for i in range(3))


def verify(root):
    light = run_path(str(root/"ch39-local-yukawa-plucker.py"))
    shell = light["shell"]
    ctx.prec = 320
    h = acb.zeta_zero(2_306_888_100_085).imag
    delta = arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
    center = 15-delta/3
    frames = {q:light["light_basis"](q,h,center) for q in (1,*CERTIFICATES)}
    for q,basis in frames.items():
        gram = basis.conjugate().transpose()*basis
        assert all((gram[i,j]-int(i==j)).contains(0)
                   and abs(gram[i,j]-int(i==j))<arb("1e-35") for i in range(3) for j in range(3))
        lap = shell["laplacian"](q,h)
        masses = [abs(t.real*(t.real-center)) for t in lap.eig(multiple=True)]
        assert sum(t<9 for t in masses)==3 and all(t<9 or t>9+(1-delta)**2 for t in masses)
    assert sorted(i for group in GROUPS for i in group)==list(range(19))
    u = frames[1]
    for q,(values,lower,upper) in CERTIFICATES.items():
        w=[None]*19
        for group,value in zip(GROUPS,values,strict=True):
            for i in group:
                w[i]=arb(value)
        assert min(w)>0
        v=frames[q]
        bound=1-compression_upper(u,w)*compression_upper(v,[1/a for a in w])
        assert bound>arb(lower)
        norms = [[sum(abs(u[i,j])*abs(v[i,k]) for i in range(19)) for k in range(3)] for j in range(3)]
        j,k=max(((j,k) for j in range(3) for k in range(3)), key=lambda pair:float(norms[pair[0]][pair[1]]))
        example=1-norms[j][k]**2
        assert example<arb(upper)
        print("目标荷、最小逐态泄漏的下界、基向量逐站对齐的上界、见证列（从零数）",
              q,bound.lower().str(14),example.upper().str(14),j,k,flush=True)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
