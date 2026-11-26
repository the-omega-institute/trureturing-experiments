"""共线引力波对独立标记的检验与跨面三标记实现的曲率迹三次项。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
import sympy as sp


ROOT = Path(__file__).parent
model = run_path(str(ROOT/"ch39-hinge-curvature-inverse.py"))


def plane_curvature(k, h):
    planes, pairs = model["PLANES"], model["PAIRS"]
    _,free = model["curvature_basis"]()
    entries = []
    for a,b in pairs:
        i,j = planes[a]
        c,d = planes[b]
        entries.append(sp.expand((k[i]*k[c]*h[j,d]+k[j]*k[d]*h[i,c]
                                 - k[i]*k[d]*h[j,c]-k[j]*k[c]*h[i,d])/2))
    return sp.Matrix([entries[i] for i in free])


def diagonal_test(matrix, sines, areas):
    # 每个无序读数代表两种反向铰链。先使五个非度量方向满足线性场方程。
    k = sp.Matrix([1,1,-1,-1])+sp.I*sp.ones(4,1)/sp.sqrt(5)
    e,f = sp.Matrix([1,-1,0,0]),sp.Matrix([0,0,1,-1])
    waves = [plane_curvature(k,e*f.T+f*e.T),plane_curvature(k,e*e.T-f*f.T)]
    deweight = sp.diag(*(1/v for v in sines))
    heavy = deweight*matrix[:,20:]
    bare = deweight*matrix[:,:20]*sp.Matrix.hstack(*waves)
    corrected = (bare-heavy*(heavy.T*heavy).inv()*heavy.T*bare).applyfunc(sp.simplify)
    assert (heavy.T*corrected).applyfunc(sp.simplify) == sp.zeros(5,2)
    assert (areas.T*corrected).applyfunc(sp.simplify) == sp.zeros(15,2)
    values = sp.zeros(2,2)
    for i,(cross,plus) in enumerate(corrected.tolist()):
        col = 0 if sines[i]==1 else 1
        # 第三腿动量 −2k 贡献四倍，成对反向读数再贡献二倍。
        values[0,col] += 8*cross**2*plus
        values[1,col] += 8*plus**3
    det = sp.simplify(values.det())
    assert det != 0
    print("共线三波的两类角亏立方矩阵",values.applyfunc(sp.simplify))
    print("两个在壳必要条件的行列式",det)


def complex_triangle():
    n,m = sp.Matrix([1,sp.I,0,0]),sp.Matrix([0,0,1,sp.I])
    momenta = [n,m,-n-m]
    polarizations = [m.conjugate(),n.conjugate(),n.conjugate()-m.conjugate()]
    curvature = []
    for k,e in zip(momenta,polarizations):
        assert (k.T*k)[0]==0 and (e.T*e)[0]==0 and (k.T*e)[0]==0
        bivector = sp.Matrix([k[i]*e[j]-k[j]*e[i] for i,j in model["PLANES"]])
        curvature.append(bivector*bivector.T)
    result = 6*sp.trace(curvature[0]*curvature[1]*curvature[2])
    assert result == 384
    return result


def verify():
    matrix,sines,areas = model["readouts"]()
    diagonal_test(matrix,sines,areas)
    operators = model["operator_coefficients"](matrix)
    bound = sum((v.applyfunc(abs) for v in operators),sp.zeros(6))
    assert sp.trace(bound**3) == sp.Rational(72152,3)
    response = complex_triangle()
    shell = run_path(str(ROOT/"ch39-arc-relative-determinant.py"))
    ctx.prec=192
    eta = [shell["phases"](acb.zeta_zero(n).imag)[724].real for n in (11366,11431)]
    b = 2*arb(3).sqrt()-3
    a = 1-b
    epsilon = arb(3)/72155
    difference = -epsilon*(eta[0]-eta[1])*(a*b)**3*int(response)/(3*(8*arb.pi())**3)
    assert -arb("4e-9") < difference < -arb("3e-9")
    print("给定复在壳三波与概率归一，两个实际零点的协变六导数顶点差",difference.str(14))


if __name__=="__main__":
    verify()
