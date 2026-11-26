"""实际零点之间随轻谱搬运三重性的有限链及两只钟的轻子空间曲率。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


study = run_path(str(Path(__file__).with_name("ch39-three-leg-yukawa.py")))
adjoint = study["adjoint"]
shell = study["shell"]
charges = (1,-4,2,-3,6)


def shell_data(height, center, charge):
    lap = shell["laplacian"](charge,height)
    eigen, vectors = lap.eig(right=True,algorithm="rump")
    basis = acb_mat(19,19)
    for j in range(19):
        assert all(not eigen[j].overlaps(r) for r in eigen[:j])
        size = sum(abs(vectors[i,j])**2 for i in range(19)).sqrt()
        for i in range(19):
            basis[i,j] = vectors[i,j]/size
    low = [j for j,r in enumerate(eigen) if abs(r.real*(r.real-center))<9]
    high = [j for j,r in enumerate(eigen) if abs(r.real*(r.real-center))>9]
    assert len(low)==3 and len(high)==16
    derivatives = []
    for chord in (724,804):
        tangent = acb_mat(19,19)
        for i,j,d,_ in shell["geo"]["edges"]:
            if d == chord:
                tangent[i,j] = acb(0,charge)*lap[i,j]
                tangent[j,i] = tangent[i,j].conjugate()
        derivatives.append(adjoint(basis)*tangent*basis)
    def gram(u,v):
        return sum((u[i,j].conjugate()*v[i,j]).real/(eigen[j].real-eigen[i].real)**2
                   for i in high for j in low)
    a,b = derivatives
    metric = acb_mat([[gram(x,y) for y in derivatives] for x in derivatives])
    curvatures = [2*sum((a[i,j].conjugate()*b[i,j]).imag/(eigen[j].real-eigen[i].real)**2
                        for i in high) for j in low]
    modes = study["edge"]["mode_columns"](study["edge"]["filters"]["spectra"](charge,height,center))
    return modes,metric,curvatures


def verify():
    ctx.prec = 448
    frames = []
    for n in (22027,28325):
        height = acb.zeta_zero(n).imag
        neutral = shell["laplacian"](0,height)
        center = sum(neutral[i,i].real for i in range(19))/18
        metrics,curvatures,frame = [],[],[]
        for q in charges:
            basis,g,f = shell_data(height,center,q)
            frame.append(basis)
            metrics.append(g)
            curvatures.extend(f)
        metric = sum(metrics,acb_mat(2,2))
        least = min(v.real.lower() for v in metric.eig(algorithm="rump"))
        print(n,"五荷共同投影变化的二阶型最小本征值、贝里曲率范围",
              least.str(9), min(z.lower() for z in curvatures).str(9),
              max(z.upper() for z in curvatures).str(9))
        assert least > arb(".0617" if n==22027 else ".1160")
        frames.append(frame)
    for q,source,target in zip(charges,*frames):
        overlap = adjoint(target)*source
        assert all(overlap[i,j].contains(0) for i in range(3) for j in range(3) if i!=j)
        defects = [1-abs(overlap[i,i])**2 for i in range(3)]
        assert all(0<d<1 for d in defects)
        print(q,"三个主角的正弦平方",*(v.str(8) for v in defects))


def triangle_phases(nodes=(11661,11934,12207,22027,28325), *, pairs=21, winner=(11934,12207,28325)):
    print("允许的零点基点", nodes)
    frames = {}
    for n in nodes:
        h = acb.zeta_zero(n).imag
        neutral = shell["laplacian"](0,h)
        c = sum(neutral[i,i].real for i in range(19))/18
        frames[n] = {q:study["edge"]["mode_columns"](study["edge"]["filters"]["spectra"](q,h,c))
                     for q in (1,-4,2)}
    best = {q:arb(0) for q in (1,-4,2)}
    winners = {}
    best_phases = {}
    triangles = list(combinations(nodes,3))
    for triangle in triangles:
        for q in best:
            a,b,c = [frames[n][q] for n in triangle]
            path = adjoint(a)*c*adjoint(c)*b*adjoint(b)*a
            phases = [path[k,k].arg() for k in range(3)]
            assert all(abs(phi)<arb(".01") for phi in phases)
            width = max(phi.upper() for phi in phases)-min(phi.lower() for phi in phases)
            eta = 2*(width/4).sin().upper()
            if eta > best[q]:
                best[q] = eta
                winners[q] = triangle
                best_phases[q] = phases
    mods = list(map(arb,(".23",".044",".0039","1")))
    coefficient = sum(a*b*c for a,b,c in combinations(mods,3))
    for q in best:
        cost = arb(".000030")/(best[q]*coefficient)
        print(q, len(triangles), "种三角形的最大去公共相位距离、取得该界的编号、必要圈数下界",
              best[q].str(10),winners[q],cost.str(9))

    s12,s23,s13 = [arb(x).sin() for x in (".22704",".04275",".00377")]
    c12,c23,c13 = [(1-s*s).sqrt() for s in (s12,s23,s13)]
    delta = arb("1.139")
    phase = acb(0,delta).exp()
    jarlskog = s12*s23*s13*c12*c23*c13**2*delta.sin()
    cs_mod = abs(c12*c23-s12*s23*s13*phase)
    coefficient = sum(a*b*c for a,b,c in combinations((s12*c13,s23*c13,s13,cs_mod),3))
    if pairs == 21:
        assert jarlskog > 37*max(best.values())*coefficient
    # 同一行、同一列各自只有一份虚部平方范数预算。
    budget = 2*pairs-1
    eta = budget*max(best.values())
    quartet_angle = (phase*(c12*c23-s12*s23*s13*phase).conjugate()).arg()
    assert arb("1.13") < quartet_angle < arb.pi()/2
    caps = []
    for j in range(128):
        lo,hi = arb(j)/128,arb(j+1)/128
        rest = eta*(1-lo**2).sqrt()
        caps.append((eta*hi/s13).asin()+(rest/(s12*c13)).asin()
                    +(rest/(s23*c13)).asin()+(eta/cs_mod).asin())
    limit = max(v.upper() for v in caps)
    assert limit < arb("1.13")
    print(budget, "圈的重定相不变量相角上界和参考相角",limit.str(10),quartet_angle.str(10))
    angles = sorted(best_phases[-4],key=float)
    width = angles[-1]-angles[0]
    assert winners[-4] == winner
    # 选相位相距最远的两条线，第三条在正反回路中抵消。
    dmat = acb_mat([[acb(0,angles[-1]).exp(),0],[0,acb(0,angles[0]).exp()]])
    target = acb_mat([[c13,s13*phase.conjugate()],[-s13*phase,c13]])
    a,d,b = cartan_angles(target)
    assert (pairs-1)*width < d < pairs*width
    print("构造三角的相位宽度", width.str(12))
    beta = ((d/pairs).sin()/abs(width.sin())).asin()
    block = dmat*rotation(beta)*adjoint(dmat)
    left,step,right = cartan_angles(block)
    assert (step-d/pairs).contains(0)
    turn = rotation(-left)*block*rotation(-right)
    reconstructed = rotation(a)*turn**pairs*rotation(b)
    assert all(abs(x)<arb("1e-35") for x in (target-reconstructed).entries())
    print("参考 CKM 的最小圈数", 2*pairs, "正反各", pairs, "次。单次中间实转角", beta.str(10))


def rotation(theta):
    s,c = theta.sin_cos()
    return acb_mat([[c,s],[-s,c]])


def cartan_angles(matrix):
    x,y = matrix[0,0],matrix[0,1]
    total = acb(x.real,y.real).arg()
    diff = acb(x.imag,-y.imag).arg()
    phase = abs(acb(x.imag,y.imag)).asin()
    return (total+diff)/2,phase,(total-diff)/2


if __name__ == "__main__":
    verify()
    triangle_phases()
    mirror = (28325,48067,58214)
    triangle_phases(mirror, pairs=14, winner=mirror)
