"""两个实际零点上，单一局域顶角要达到夸克总质量比所需的重 Higgs 份额。"""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, ctx

local = run_path(str(Path(__file__).with_name("ch39-local-light-higgs-hierarchy.py")))
adj = local["adjoint"]


def operator(left,right):
    return acb_mat([[left[k,i]*right[k,j] for k in range(19)] for i in range(3) for j in range(3)])


def fixed(x):
    return arb(float(x))


def vector(x):
    return acb_mat([[acb(fixed(z.real),fixed(z.imag))] for z in x])


def numeric(M):
    return np.array([[complex(M[i,j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])


def quadratic(x,M):
    return (adj(x)*M*x)[0,0].real


def positive_definite(M):
    # 对厄米矩阵作无主元 LDL* 消元，每个实枢轴都严格为正才接受。
    A=M+acb_mat(M.nrows(),M.ncols())
    for i in range(M.nrows()):
        pivot=A[i,i].real
        assert pivot>0
        for j in range(i+1,M.nrows()):
            for k in range(i+1,j+1):
                A[j,k] -= A[j,i]*A[k,i].conjugate()/pivot


def certificate(P,B,eye,A,D,low_target,high_target):
    p,b=numeric(P),numeric(B)
    def optimizer_at(lam):
        return np.linalg.eigh(p-lam*b)[1][:,0]
    def constraint(lam):
        z=optimizer_at(lam)
        return np.vdot(z,b@z).real
    upper=1.0
    while constraint(upper)<0:
        upper*=2
        assert upper<1e7
    lower_lam=0.
    for _ in range(55):
        mid=(lower_lam+upper)/2
        if constraint(mid)>0:
            upper=mid
        else:
            lower_lam=mid
    lam=(lower_lam+upper)/2
    lower=fixed(np.linalg.eigvalsh(p-lam*b)[0]-1e-7)
    positive_definite(P-fixed(lam)*B-lower*eye)
    # 两侧最低向量在退化点也可混合。搜索只建议向量，最终用区间形式验证可行性。
    def proposal(scale):
        return np.linalg.eigh(p-(lam*scale)*b)[1][:,0]
    x,y=proposal(.9999),proposal(1.0001)
    phase=np.vdot(x,y)
    if abs(phase)>1e-15:
        y=y*phase.conjugate()/abs(phase)
    lo,hi=0.,1.
    for _ in range(60):
        mid=(lo+hi)/2
        trial=(1-mid)*x+mid*y
        if np.vdot(trial,b@trial).real>0:
            hi=mid
        else:
            lo=mid
    candidates=[]
    for offset in (1e-4,.01,.1,1.):
        t=min(1.,hi+offset)
        z=vector((1-t)*x+t*y)
        if quadratic(z,B)>0:
            norm_down=sum(abs(entry)**2 for entry in D*z)
            norm_up=sum(abs(entry)**2 for entry in A*z)
            if norm_down>0 and low_target < (norm_up/norm_down).sqrt() < high_target:
                candidates.append((quadratic(z,P)/quadratic(z,eye)).upper())
    assert candidates
    ceiling=min(candidates)
    assert 0<ceiling-lower<arb(".0001")
    return lower,ceiling


def verify():
    ctx.prec=384
    up_values=list(map(arb,(".00000615",".00311",".8616")))
    dn_values=list(map(arb,(".0000135",".000268",".01401")))
    ratio=(sum(y*y for y in up_values)/sum(y*y for y in dn_values)).sqrt()
    low=ratio*(arb(".999")/arb("1.001"))
    eye=acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    print("中心范数比与参考窗下端",ratio.str(9),low.str(9))
    for number in (22027,28325):
        height=acb.zeta_zero(number).imag
        L0=local["local"]["shell"]["laplacian"](0,height)
        c=sum(L0[i,i].real for i in range(19))/18
        V={q:local["basis"](q,height,c) for q in (1,-4,2,3)}
        A=operator(V[1],V[-4])
        D=operator(V[1],V[2]).conjugate()
        P=eye-V[3]*adj(V[3])
        B=adj(A)*A/low**2-adj(D)*D
        lo,hi=certificate(P,B,eye,A,D,low,ratio*(arb("1.001")/arb(".999")))
        assert lo>arb(".99858")
        print(number,"达到窗下端所需重 Higgs 范数平方份额的最小值",lo.str(10),hi.str(10))


if __name__=="__main__":
    verify()
