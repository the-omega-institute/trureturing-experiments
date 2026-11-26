"""把不同超荷的原壳算子压缩到三个原始 Gram 源，检验共同味标架与 CP。"""
from pathlib import Path
from runpy import run_path
from itertools import combinations

from flint import acb, acb_mat, arb, ctx


def identity(n):
    return acb_mat([[int(i==j) for j in range(n)] for i in range(n)])


def adj(matrix):
    return matrix.conjugate().transpose()


def family_space(gram):
    # 求 g1,g2,g3 在四点空间内的共同实法向量。沿首分量选正号。
    normal = [arb(1)]
    for i in (1,2,3):
        normal.append(-sum(gram[i][k]*normal[k] for k in range(i))/gram[i][i])
    length = sum(x*x for x in normal).sqrt()
    n = [x/length for x in normal]
    v = acb_mat([[x+int(j==0)] for j,x in enumerate(n)])
    reflection = identity(4)-v*adj(v)*(1/(1+n[0]))
    frame = acb_mat([[reflection[i,j] for j in (1,2,3)] for i in range(4)])
    assert all((a-b).contains(0) for a,b in zip((adj(frame)*frame).entries(),identity(3).entries()))
    overlap = n[0]**2+n[1]*n[2]+n[2]*n[3]+n[3]*n[1]
    assert abs(overlap)<1
    return frame, overlap


def projections(matrix):
    levels = sorted((v.real for v in matrix.eig(algorithm="rump")),key=lambda v:float(v.mid()))
    assert levels[0]>0 and levels[0]<levels[1]<levels[2]
    one = identity(3)
    projections = []
    for i,e in enumerate(levels):
        p = one
        for j,f in enumerate(levels):
            if j!=i:
                p = (p*(matrix-one*f))*(1/(e-f))
        projections.append(p)
    return levels,projections


def mixing(up,down):
    up_e, up_p = projections(up)
    down_e, down_p = projections(down)
    probabilities = [[(p*q).trace().real for q in down_p] for p in up_p]
    assert all(0<t<1 for row in probabilities for t in row)
    # Im(V11 V22 conj(V12 V21)) 与交换子行列式的绝对值独立求得。
    j = (up_p[0]*down_p[0]*up_p[1]*down_p[1]).trace().imag
    divisor = arb(2)
    for levels in (up_e,down_e):
        for a,b in combinations(levels,2):
            divisor *= b-a
    det_j = abs((up*down-down*up).det().imag)/divisor
    assert abs(j).overlaps(det_j)
    s13 = probabilities[0][2].sqrt()
    c13 = (1-s13*s13).sqrt()
    min_mod = min(p.sqrt().lower() for row in probabilities for p in row)
    return [probabilities[0][1].sqrt()/c13,probabilities[1][2].sqrt()/c13,s13,j], (up_e,down_e), min_mod


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    ckm = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ctx.prec = 320
    gram = native["gram_cholesky"](root)
    frame, overlap = family_space(gram)
    positions = [shell["soft"]["pole"],*shell["boundary_sites"]()]
    assert positions==[14,8,7,6]
    height = acb.zeta_zero(2306888100085).imag
    print("三源法向量与其三循环像的内积",overlap.str(14))
    direct, square, zero_energy = [], [], []
    labels = shell["geo"]["sites"]
    rotate = [labels.index((y,z,x)) for x,y,z in labels]
    assert rotate[positions[0]] == positions[0] and set(rotate[i] for i in positions[1:])==set(positions[1:])
    for charge in (-4,2):
        matrix = shell["laplacian"](charge,height)
        assert all((matrix[i,j]-matrix[rotate[i],rotate[j]]).contains(0) for i in range(19) for j in range(19))
        assert all(matrix[positions[0],s].imag.is_zero() and matrix[positions[0],s].real<0 for s in positions[1:])
        for operator, target in ((matrix,direct),(matrix*matrix,square),(matrix.inv(),zero_energy)):
            block = acb_mat([[operator[i,j] for j in positions] for i in positions])
            target.append(adj(frame)*block*frame)
        diag = (matrix*matrix)[positions[0],positions[1]].imag
        print("荷、二步极点至首站虚部",charge,diag.str(11))
        assert not diag.contains(0)
    ref = [v+arb(0,ckm["sine_radius"].upper()) for v in ckm["centers"][:3]]
    for name,(up,down) in zip(("一次","二次","零能逆算子"),(direct,square,zero_energy)):
        results, (up_e,down_e), min_mod = mixing(up,down)
        print(name,"九个模中的最小下界", min_mod.str(12))
        assert min_mod>arb(".0138")>ref[2]
        assert all(not p.overlaps(r) for p,r in zip(results[:3],ref,strict=True))
        if name=="一次":
            assert results[-1].contains(0)
        elif name=="二次":
            assert arb(".00030018913449")<results[-1]<arb(".00030018913450")
        else:
            assert arb("-.000005553")<results[-1]<arb("-.000005552")
        print(name,"上、下三阶压缩谱",*[x.str(10) for x in up_e],*[x.str(10) for x in down_e])
        print(name,"s12、s23、s13、有符号 J",*[v.str(14) for v in results],flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
