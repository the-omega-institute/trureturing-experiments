"""五次全局反例，以及冻结高度的六次谱作用局部校准。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx


def chebyshev(x):
    u, v = (arb(1),arb(0),arb(0)), (x,arb(1),arb(0))
    out = [v]
    for _ in range(5):
        cur = (2*x*v[0]-u[0], 2*v[0]+2*x*v[1]-u[1], 4*v[1]+2*x*v[2]-u[2])
        u,v = v,cur
        out.append(cur)
    return out


def verify(replay):
    ctx.prec=256
    s=run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump=run_path(str(replay/"ch39-native-light-projector-pump.py"))
    adj=pump["adjoint"]
    bases=list(pump["rotation_bases"](s))
    def compute(q,n,with_jet):
        h=acb.zeta_zero(n).imag
        lap=s["laplacian"](q,h)
        result=[arb(0) for _ in range(6)]
        if with_jet:
            # 此处坐标是作用在荷 q 后的两条有向钟角，每个相位因子 e^{i phi_l}。
            _,terms=pump["clock_matrices"](s,q)
            amplitudes=[m*s["geo"]["phases"](h)[d]**q for m,d in zip(terms,(724,804))]
            first=[acb(0,1)*(m-adj(m)) for m in amplitudes]
            second=[-(m+adj(m)) for m in amplitudes]
            grad=[[arb(0) for _ in range(6)] for k in range(2)]
            hess=[[arb(0) for _ in range(6)] for k in range(3)]
        for basis in bases:
            H=adj(basis)*lap*basis
            if with_jet:
                levels,vec=H.eig(right=True,algorithm="rump")
            else:
                levels=H.eig(algorithm="rump")
            order=sorted(range(len(levels)),key=lambda i:float(levels[i].real.mid()))
            assert all(v.is_finite() and v.imag.contains(0) for v in levels)
            lam=[levels[i].real for i in order]
            assert all(a<b for a,b in zip(lam,lam[1:]))
            jets=chebyshev(2*lam[0]/9-1)
            for k,(value,fp,fpp) in enumerate(jets):
                result[k]+=value
            if not with_jet:
                continue
            vectors=[]
            for j in order:
                v=acb_mat([[vec[i,j]] for i in range(H.nrows())])
                vectors.append(v/sum(abs(t)**2 for t in v.entries()).sqrt())
            ground=vectors[0]
            deriv=[adj(basis)*x*basis for x in first]
            h2=[adj(basis)*x*basis for x in second]
            couplings=[[(adj(v)*D*ground)[0,0] for v in vectors] for D in deriv]
            dx=[2*a[0].real/9 for a in couplings]
            dd=[]
            for i,j in ((0,0),(0,1),(1,1)):
                direct=(adj(ground)*h2[i]*ground)[0,0].real if i==j else arb(0)
                mixed=sum(2*(couplings[i][k].conjugate()*couplings[j][k]).real/(lam[0]-lam[k]) for k in range(1,len(lam)))
                dd.append(2*(direct+mixed)/9)
            for k,(_,fp,fpp) in enumerate(jets):
                for i in range(2):
                    grad[i][k]+=fp*dx[i]
                for m,(i,j) in enumerate(((0,0),(0,1),(1,1))):
                    hess[m][k]+=fp*dd[m]+fpp*dx[i]*dx[j]
        return (result,grad,hess) if with_jet else result
    cols=[];gradrows=[];curv=[]
    for q,n,others in ((-4,2345666,(31681,1162532)),(2,964890,(239010,500730,2634783,23551))):
        values,grad,hess=compute(q,n,True)
        gradrows+=grad;curv.append(hess)
        for other in others:
            new=compute(q,other,False)
            cols.append([new[j]-values[j] for j in range(5)])
    A=arb_mat([[v[i] for v in cols] for i in range(5)]+[[1]*6])
    R=A.inv()
    bary=[R[i,5] for i in range(6)]
    distances=[bary[i]/sum(R[i,j]**2 for j in range(5)).sqrt() for i in range(6)]
    assert all(v>0 for v in bary)
    print("五维单形增广行列式",A.det().str(13))
    print("重心坐标",[v.str(15) for v in bary])
    print("到各面的距离",[v.str(15) for v in distances])

    assert abs(A.det()) > arb("1.88198e-13")
    assert all(d > arb("5.37486226e-7") for d in distances)
    assert all(distances[2] < d for i,d in enumerate(distances) if i!=2)
    gradient = arb_mat(gradrows)
    def sub(columns):
        return arb_mat([[row[i] for i in columns] for row in gradrows])
    # 四次及以下的四个非恒定系数无非零驻点解。
    quartic_det = sub(range(4)).det()
    assert abs(quartic_det) > 0
    print("四次公共驻点方程的行列式", quartic_det.str(14))
    # 分别指定 (c1,c2)=(-1,0) 与 (c1,c6)=(-1,0)，其余由四个驻点条件唯一解出。
    for degree, unknown in ((6, (2,3,4,5)), (5, (1,2,3,4))):
        coefficient = arb_mat(6,1)
        coefficient[0,0] = -1
        mat = sub(unknown)
        assert abs(mat.det()) > 0
        solved = mat.solve(sub((0,)))
        for j,index in enumerate(unknown):
            coefficient[index,0] = solved[j,0]
        assert all(v.contains(0) for v in (gradient*coefficient).entries())
        assert abs(coefficient[degree-1,0]) > 0
        print("次数、系数", degree,[v.str(15) for v in coefficient.entries()])
        for k,rows in enumerate(curv):
            a,b,d=(arb_mat(rows)*coefficient).entries()
            least=(a+d-((a-d)**2+4*b*b).sqrt())/2
            assert least > arb(".013")
            print("部门、Hessian 最小本征值",(-4,2)[k],least.str(16))


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
