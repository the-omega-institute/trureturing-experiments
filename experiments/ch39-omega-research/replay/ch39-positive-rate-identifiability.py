"""检验实际零点平方读数的连续前缀。用全正概率碰撞给下界，用二阶子式及正定性给上界。"""
import os
from pathlib import Path
from runpy import run_path
import sys

# 稠密小矩阵的浮点搜索可在单线程重现。严格结论来自最后的 Arb 正定检验。
for var in ("OPENBLAS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","MKL_NUM_THREADS"):
    os.environ[var]="1"
import numpy as np
from flint import arb, arb_mat, ctx


COUNT = 260
ZERO_ROWS = (1,3,8,11,15,17,18,20,21,23,25,26,27,29,30,33,36,37,40,43,44,45,48,49,50,53,56,58)
SMALL = (2,3,29)


def collision(rows):
    fixed = set(SMALL)
    free = [j for j in range(32) if j not in fixed]
    restrictions = [[arb(1)]*32]+[rows[i-1] for i in ZERO_ROWS]
    epsilon = arb(1)/10**5
    mat = arb_mat([[row[j] for j in free] for row in restrictions])
    rhs = arb_mat([[arb(k==0)-epsilon*sum(row[j] for j in fixed)] for k,row in enumerate(restrictions)])
    assert not mat.det().contains(0)
    values = mat.solve(rhs).entries()
    assert all(v>epsilon for v in values) and len(free)==len(ZERO_ROWS)+1
    print("前缀碰撞、共同消零行数、自由概率的最小值",58,len(ZERO_ROWS),min(z.lower() for z in values).str(15),flush=True)
    assert 58-len(ZERO_ROWS)==30


def partial_transpose(matrix):
    d=32*31
    return matrix.reshape(32,31,32,31).transpose(0,3,2,1).reshape(d,d)


def change_of_basis(rows):
    # 和为零切空间使用 e_j-e_31 为有理基。两个右乘的可逆矩阵专为数值条件而选。
    m = np.asarray(rows,dtype=np.float64)
    raw = m[:,:31]-m[:,31,None]
    _,upper = np.linalg.qr(m,mode="reduced")
    _,tangent = np.linalg.qr(raw,mode="reduced")
    X,Y = np.linalg.inv(upper),np.linalg.inv(tangent)
    a,b = m@X,raw@Y
    D=(a[:,:,None]*b[:,None,:]).reshape(len(rows),-1)
    return X,Y,D


def propose(D):
    # 每次正特征方向都是当前投影的修正量，只在浮点中选证书。
    n,d=D.shape
    Q=np.linalg.qr(D.T,mode="complete")[0]
    U,K=Q[:,:n],Q[:,n:]
    size=K.shape[1]
    restricted=np.zeros((size,size))
    coefficient=np.zeros_like(restricted)
    prev_h,prev_c=restricted,coefficient
    t=1
    for turn in range(900):
        ev,v=np.linalg.eigh(restricted)
        if ev[0]>.025:
            break
        active=ev<1
        H=(v[:,active]*(1-ev[active]))@v[:,active].T
        lift=K@H@K.T
        next_h=restricted+(H-K.T@partial_transpose(lift)@K)/2
        next_c=coefficient+H/2
        t2=(1+np.sqrt(1+4*t*t))/2
        beta=(t-1)/t2
        restricted=next_h+beta*(next_h-prev_h)
        coefficient=next_c+beta*(next_c-prev_c)
        prev_h,prev_c,t=next_h,next_c,t2
    else:
        raise ArithmeticError("浮点搜索在预设迭代内未找到正定证书")
    M=K@coefficient@K.T
    M=(M+M.T)/2
    F=M-partial_transpose(M)
    P=U@U.T
    C=np.linalg.lstsq(D.T,-F+(F@P+P)/2,rcond=None)[0]
    H=F+D.T@C+C.T@D
    L=np.linalg.cholesky((H+H.T)/2)
    N=np.linalg.inv(L)
    print("待验前缀、张量维数、浮点提议迭代数",n,d,turn,flush=True)
    return M,C,N


def ball_matrix(matrix):
    return arb_mat(matrix.tolist())


def upper(rows):
    X,Y,D=change_of_basis(rows)
    M,C,N=propose(D)
    # NumPy 只筛选矩阵。将浮点数当成精确 dyadic 有理数，用严格的 R 而非浮点 D 验证。
    ctx.prec=96
    X,Y=ball_matrix(X),ball_matrix(Y)
    assert not X.det().contains(0) and not Y.det().contains(0)
    a=arb_mat(rows)*X
    b=arb_mat([[r[j]-r[-1] for j in range(31)] for r in rows])*Y
    exact=arb_mat([[a[i,j]*b[i,k] for j in range(32) for k in range(31)] for i in range(len(rows))])
    # 分开转为精确球，再相减，令秩一向量上恒零的二阶子式不受浮点减法影响。
    F=ball_matrix(M)-ball_matrix(partial_transpose(M))
    W=exact.transpose()*ball_matrix(C)
    cert=F+W+W.transpose()
    pre=ball_matrix(N)
    checked=pre*cert*pre.transpose()
    d=32*31
    deviation=max(sum(abs(checked[i,j]-int(i==j)) for j in range(d)).upper() for i in range(d))
    print("整张严格矩阵预乘后的最大绝对行和误差",deviation.str(13),flush=True)
    assert deviation<arb(".01")


def verify(root):
    source=run_path(str(root/"ch39-proton-detection-prefix.py"))
    _,_,rows=source["samples"](root,COUNT)
    collision(rows)
    upper(rows)


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
