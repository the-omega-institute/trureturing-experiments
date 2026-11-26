"""十九点壳上的逐站点折叠，三个轻 Higgs 任意叠加的味谱界。"""
from itertools import combinations, combinations_with_replacement
from pathlib import Path
from runpy import run_path

from flint import acb_mat, acb, arb, ctx


local = run_path(str(Path(__file__).with_name("ch39-local-yukawa-plucker.py")))
basis = local["light_basis"]


def adjoint(matrix):
    return matrix.conjugate().transpose()


def tensors(left, right, higgs):
    blocks = []
    for k in range(3):
        local_right = acb_mat([[higgs[i,k]*right[i,j] for j in range(3)] for i in range(19)])
        blocks.append(left.transpose()*local_right)
    return blocks


def gram(blocks):
    return acb_mat([[sum(a[i,j].conjugate()*b[i,j] for i in range(3) for j in range(3))
                     for b in blocks] for a in blocks])


def wedge_bounds(blocks, metric):
    pairs = list(combinations(range(3),2))
    quadratic = list(combinations_with_replacement(range(3),2))
    rows = []
    for i,j in pairs:
        for k,l in pairs:
            def coeff(a,b):
                A,B = blocks[a],blocks[b]
                return A[i,k]*B[j,l]-A[i,l]*B[j,k]
            rows.append([coeff(a,b) if a==b else (coeff(a,b)+coeff(b,a))/arb(2).sqrt()
                         for a,b in quadratic])
    W = acb_mat(rows)
    eigen = (adjoint(W)*W).eig(algorithm="rump")
    minimum = min(z.real.lower() for z in eigen)
    maximum = max(z.real.upper() for z in metric.eig(algorithm="rump"))
    assert minimum>0 and maximum>0
    return minimum.sqrt()/(arb(3).sqrt()*maximum)


def extreme_ratio(G,H):
    assert all(z.real>0 for z in H.eig(algorithm="rump"))
    eig = (H.inv()*G).eig(algorithm="rump")
    return min(z.real.lower() for z in eig).sqrt(),max(z.real.upper() for z in eig).sqrt()


def verify():
    ctx.prec=384
    for number in (22027,28325):
        h=acb.zeta_zero(number).imag
        L0=local["shell"]["laplacian"](0,h)
        center=sum(L0[i,i].real for i in range(19))/18
        V={q:basis(q,h,center) for q in (1,-4,2,-3,6,3)}
        up=tensors(V[1],V[-4],V[3])
        down=tensors(V[1],V[2],V[3].conjugate())
        leptons=tensors(V[-3],V[6],V[3].conjugate())
        G,H,K=[gram(blocks) for blocks in (up,down,leptons)]
        print(number, "上/下、下/带电轻子的 Frobenius 范数比区间",
              *[x.str(9) for x in (*extreme_ratio(G,H.conjugate()),*extreme_ratio(H,K))])
        floors=[wedge_bounds(blocks,metric) for blocks,metric in ((up,G),(down,H),(leptons,K))]
        print("任意三个轻 Higgs 混合的次大/最大奇异值下界",*[v.str(9) for v in floors])


if __name__=="__main__":
    verify()
