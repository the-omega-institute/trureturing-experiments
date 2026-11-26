"""联合近真空零点的轮换特征标和一跳协变整代传递。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


def adjoint(matrix):
    return matrix.conjugate().transpose()


def generator(shell):
    """以裸站点轮换实现 C3，并对两只钟的有向边和软辐边逐条核对。"""
    geo=shell["geo"]
    sites=geo["sites"]
    perm=[sites.index((b,c,a)) for a,b,c in sites]
    assert all(perm[perm[perm[i]]]==i for i in range(19))
    pole=shell["soft"]["pole"]
    assert sites[pole]==(17,17,17)
    orbit=shell["soft"]["star"]["orbit"]
    assert perm[pole]==pole and set(perm[i] for i in orbit)==set(orbit)
    edges={(i,j):length for i,j,length,_ in geo["edges"]}
    for (i,j),length in edges.items():
        ri,rj=perm[i],perm[j]
        assert edges[tuple(sorted((ri,rj)))]==length
        if length in (724,804):
            assert ri<rj
        # 原生相位的实部由内积给出，虚部由以极点为首行的有向行列式给出。
        assert geo["dot"](sites[ri],sites[rj])==geo["dot"](sites[i],sites[j])
        assert geo["det"](sites[pole],sites[ri],sites[rj])==geo["det"](sites[pole],sites[i],sites[j])
    return acb_mat([[int(i==perm[j]) for j in range(19)] for i in range(19)]),edges


def character_order(frame,S):
    omega=acb(-arb(1)/2,arb(3).sqrt()/2)
    chars=[acb(1),omega,omega.conjugate()]
    matrix=adjoint(frame)*S*frame
    order=[]
    for j in range(3):
        matches=[k for k,z in enumerate(chars) if (matrix[j,j]-z).contains(0)]
        assert len(matches)==1 and abs(matrix[j,j]-chars[matches[0]])<arb("1e-30")
        order.append(matches[0])
    return order


def verify(root):
    light=run_path(str(root/"ch39-local-yukawa-plucker.py"))
    prev=run_path(str(root/"ch39-joint-zero-diagonal-light-filters.py"))
    shell=light["shell"]
    ctx.prec=320
    height=acb.zeta_zero(2_306_888_100_085).imag
    delta=arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    center=15-delta/3
    S,edges=generator(shell)
    eye=acb_mat([[int(i==j) for j in range(19)] for i in range(19)])
    chars={}
    frames={}
    for q in (1,2,-4,-3,6):
        lap=shell["laplacian"](q,height)
        assert all(z.contains(0) and abs(z)<arb("1e-35") for z in (S*lap-lap*S).entries())
        frame=light["light_basis"](q,height,center)
        order=character_order(frame,S)
        chars[q]=order
        frames[q]=frame
        print("荷、本征线按返回顺序的轮换特征标与重数",q,order,[order.count(i) for i in range(3)],flush=True)
    blocked=[(i,j) for i in range(19) for j in range(19)
             if i!=j and tuple(sorted((i,j))) not in edges]
    assert len(blocked)==72
    u=frames[1]
    p=u*adjoint(u)
    pc=eye-p
    for q in (2,-4):
        assert sorted(chars[1])==sorted(chars[q])==[0,1,2]
        cols=[chars[q].index(k) for k in chars[1]]
        v=acb_mat([[frames[q][i,k] for k in cols] for i in range(19)])
        qc=eye-v*adjoint(v)
        prev["certify_positive"](acb_mat([[qc[i,j]*p[j,i] for j in range(19)] for i in range(19)]))
        base=v*adjoint(u)
        H=acb_mat([[qc[i,k]*pc[l,j] for k,l in blocked] for i,j in blocked])
        prev["certify_positive"](H-arb(".033")*acb_mat([[int(i==j) for j in range(72)] for i in range(72)]))
        # 将三个特征标的补块分开认证，可同时覆盖每条线的任意传输相位。
        split=H.solve(acb_mat([[-v[i,a]*u[j,a].conjugate() for a in range(3)] for i,j in blocked]))
        bounds=[]
        for a in range(3):
            sparse=acb_mat(19,19)
            for index,(i,j) in enumerate(blocked):
                sparse[i,j]=split[index,a]
            correction=qc*sparse*pc
            bounds.append(sum(abs(z)**2 for z in correction.entries()).sqrt())
        total=sum(bounds)
        assert total<(arb(".893") if q==2 else arb(".915"))
        print("三个特征标补块的 Frobenius 范数及总和",q,*[x.str(10) for x in (*bounds,total)],flush=True)
        coeff=acb_mat([[sum(split[i,j] for j in range(3))] for i in range(72)])
        L=acb_mat(19,19)
        for (i,j),z in zip(blocked,coeff.entries(),strict=True):
            L[i,j]=z
        B=qc*L*pc
        T=base+B
        # 对称平均投影到精确的交织空间。各次轮换保留原边支撑和目标基配对。
        averaged=(T+S*T*adjoint(S)+adjoint(S)*T*S)/3
        residual=averaged-base
        norm=sum(abs(z)**2 for z in residual.entries()).sqrt()
        assert norm<1
        assert all(averaged[i,j].contains(0) and abs(averaged[i,j])<arb("1e-30") for i,j in blocked)
        for diff in (S*averaged-averaged*S,averaged*u-v,adjoint(averaged)*v-u):
            assert all(z.contains(0) and abs(z)<arb("1e-30") for z in diff.entries())
        print("目标荷、按源特征标选取的目标列、协变补块的 Frobenius 范数",q,cols,norm.str(13),flush=True)


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    verify(root)
