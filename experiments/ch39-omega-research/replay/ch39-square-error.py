"""前260个平方读数的统一逆 Lipschitz 上界，以及严格正概率间的可达误差放大。"""
import os
from pathlib import Path
from runpy import run_path
import sys

for variable in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

from flint import arb, arb_mat, ctx
import numpy as np


DENOMINATOR = 10**8
FIRST = (2147676,1776354,3000100,847318,240639,100,4325440,5921407,
         100,1064575,100,1884289,8175539,6091111,1610387,2393437,
         1866608,2560617,8763097,8248908,4082027,5665004,1513015,6689142,
         3200723,6859819,4199951,6779437,100,100,72398,20482)
SECOND = (2147849,1776190,100,3401145,242799,982,4257136,5969289,
          1083406,87133,2259171,100,7735788,6513820,1609763,2392998,
          1866000,2559825,8762340,8249161,4094102,5666706,1511216,6691814,
          3200269,6859471,4200244,6777524,64883,18576,100,100)


def near_pair(rows):
    assert sum(FIRST) == sum(SECOND) == DENOMINATOR and min(*FIRST, *SECOND) > 0
    difference = [arb(a-b)/DENOMINATOR for a,b in zip(FIRST,SECOND)]
    total = [arb(a+b)/DENOMINATOR for a,b in zip(FIRST,SECOND)]
    data = [abs(sum(t*v for t,v in zip(row,difference))*sum(t*v for t,v in zip(row,total))) for row in rows]
    top = 87
    assert all(v < data[top] for i,v in enumerate(data) if i != top)
    displacement = arb(3)/100
    assert max(abs(a-b) for a,b in zip(FIRST,SECOND)) == 3_000_000
    assert data[top] < arb("2.471e-7") and displacement/data[top] > 121456
    tv = sum(abs(a-b) for a,b in zip(FIRST,SECOND))
    assert tv == 12_940_528
    print("近碰撞的最大逐项概率差、总变差、最大平方读数差及编号、放大比",
          displacement.str(10), (arb(tv)/(2*DENOMINATOR)).str(10),
          data[top].str(18), top+1, (displacement/data[top]).str(16), flush=True)


def metric_proposal(D, F):
    # 在 D 的核上 F 为正。用 Schur 补的绝对行和选择测量权，再对候选度量预条件。
    size = D.shape[0]
    q = np.linalg.qr(D.T, mode="complete")[0]
    span, null = q[:,:size], q[:,size:]
    correction = np.linalg.solve(null.T@F@null, null.T@F@span)
    inverse = np.linalg.inv(D@span)
    schur = inverse.T @ (span.T@F@(span-null@correction)) @ inverse
    weights = np.ceil(np.maximum(1, np.sum(np.abs(schur), axis=1))).astype(np.int64)
    metric = D.T@(weights[:,None]*D) + 1.42*F
    factor = np.linalg.cholesky((metric+metric.T)/2)
    return weights.tolist(), np.linalg.inv(factor)


def upper(rows, previous):
    X,Y,D = previous["change_of_basis"](rows)
    M,_,_ = previous["propose"](D)
    flip = previous["partial_transpose"](M)
    weights, pre = metric_proposal(D, M-flip)
    ctx.prec = 112
    exact_matrix = previous["ball_matrix"]
    X,Y = exact_matrix(X),exact_matrix(Y)
    assert not X.det().contains(0) and not Y.det().contains(0)
    left = arb_mat(rows)*X
    right = arb_mat([[r[j]-r[-1] for j in range(31)] for r in rows])*Y
    exact_D = arb_mat([[left[i,j]*right[i,k] for j in range(32) for k in range(31)] for i in range(len(rows))])
    weighted = arb_mat([[weights[i]*exact_D[i,j] for j in range(992)] for i in range(len(rows))])
    assert min(weights)>0
    metric = exact_D.transpose()*weighted + (arb(71)/50)*(exact_matrix(M)-exact_matrix(flip))
    pre = exact_matrix(pre)
    checked = pre*metric*pre.transpose()
    error = max(sum(abs(checked[i,j]-int(i==j)) for j in range(992)).upper() for i in range(992))
    assert error < arb(".01")
    # ||p||_1=||p'||_1=1 对概率成立；提取式只使用两者坐标和各为一。
    sums = [sum(X[i,j] for i in range(32)) for j in range(32)]
    tangent = [[Y[i,j] for j in range(31)] for i in range(31)]
    tangent.append([-sum(Y[i,j] for i in range(31)) for j in range(31)])
    extract = arb_mat([[s*t/2 for s in sums for t in row] for row in tangent])
    read = extract*pre.transpose()
    lengths = [sum(read[i,j]**2 for j in range(992)).sqrt().upper() for i in range(32)]
    bound = max(lengths)*(arb(sum(weights))/(1-error)).sqrt()
    print("严格度量预条件误差、测量权和、最大提取行范数、统一逐项逆界",
          error.str(15),sum(weights),max(lengths).str(15),bound.str(17),flush=True)
    assert bound < 43_904_358


def verify(root):
    previous = run_path(str(root/"ch39-positive-rate-identifiability.py"))
    samples = run_path(str(root/"ch39-proton-detection-prefix.py"))["samples"]
    _,_,rows = samples(root,previous["COUNT"])
    assert len(rows)==260
    near_pair(rows)
    upper(rows,previous)


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
