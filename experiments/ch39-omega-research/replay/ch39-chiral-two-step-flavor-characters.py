"""两步三代质量边界的轮换特征标与单一 Higgs 真空的味禁戒。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

ROOT = Path(__file__).parent
shell = run_path(str(ROOT / "ch39-three-ray-chiral-boundary.py"))
inertia = run_path(str(ROOT / "ch39-monotone-heavy-mass-ratio.py"))["count_above"]
window = run_path(str(ROOT / "ch39-chiral-common-gap-crossing.py"))


def identity(n):
    mat = acb_mat(n,n)
    for i in range(n):
        mat[i,i] = 1
    return mat


def verify():
    ctx.prec = 192
    sites = shell["geo"]["sites"]
    rot = [sites.index((b,c,a)) for a,b,c in sites]
    orbits = []
    for i,j in enumerate(rot):
        if j==i:
            orbits.append((i,))
        elif i==min(i,j,rot[j]):
            orbits.append((i,j,rot[j]))
    assert sorted(v for group in orbits for v in group)==list(range(19))
    root = acb(-1,arb(3).sqrt())/2
    bases = []
    for k in range(3):
        groups = [g for g in orbits if k==0 or len(g)==3]
        V = acb_mat(19,len(groups))
        for j,orbit in enumerate(groups):
            for t,i in enumerate(orbit):
                V[i,j] = root**(k*t)
        bases.append(V)
    eye = identity(19)
    for n in [*window["EXPECTED"][1:],12207]:
        height = acb.zeta_zero(n).imag
        neutral = shell["laplacian"](0,height)
        center = sum(neutral[i,i].real for i in range(19))/18
        for q in (1,-4,2,-3,6):
            L = shell["laplacian"](q,height)
            # 循环坐标置换保持边长、定向平行输运与三辐边权。
            assert all((L[rot[i],rot[j]]-L[i,j]).contains(0) for i in range(19) for j in range(19))
            D = L*(L-center*eye)
            cut = 81*eye - D*D
            counts = [inertia(V.conjugate().transpose()*cut*V,0) for V in bases]
            assert counts==[1,1,1]
    # 两个不同特征标的真空源可同时产生 I 与轮换矩阵的味张量。
    up = acb_mat([[1,0,0],[0,2,0],[0,0,4]])
    cycle = acb_mat([[0,1,0],[0,0,1],[1,0,0]])
    down = identity(3)+acb(0,1)*cycle
    a,b = (Y*Y.conjugate().transpose() for Y in (up,down))
    det = (a*b-b*a).det()
    assert det.real==0 and abs(det.imag)==1080
    print("十一枚固定分界候选，五种有符号荷的三个轮换重数均为 1,1,1")
    print("双特征标例子的味对易子行列式",det)


if __name__=="__main__":
    verify()
