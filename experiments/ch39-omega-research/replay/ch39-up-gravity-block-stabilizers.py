"""右外壳三站、六站荷六块的整个稳定子分别允许轻上夸克中心。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


BLOCKS = ((14,6,7), (14,17,8,16,4,13))


def build_inputs(root):
    lines = (root/"ch39-up-gravity-block-stabilizers.tsv").read_text().splitlines()
    rows = iter(tuple(map(int,s.split())) for s in lines if s and not s.startswith("#"))
    witnesses = []
    for group in BLOCKS:
        direct = acb_mat(19,19)
        for indices in (group, tuple(i for i in range(19) if i not in group)):
            dim = len(indices)
            z = acb_mat(dim,dim)
            for i in range(dim):
                coords = next(rows)
                assert len(coords) == 2*dim
                for j in range(dim):
                    z[i,j] = acb(coords[2*j],coords[2*j+1])
            # 精确的高斯整数列作正交归一化，得到与候选接近的精确酉矩阵。
            basis = []
            for j in range(dim):
                col = acb_mat([[z[i,j]] for i in range(dim)])
                for vector in basis:
                    col -= vector*(vector.conjugate().transpose()*col)[0,0]
                norm = (col.conjugate().transpose()*col)[0,0].real.sqrt()
                assert norm > 0
                basis.append(col*(1/norm))
            for a,i in enumerate(indices):
                for b,j in enumerate(indices):
                    direct[i,j] = basis[b][a,0]
        witnesses.append(direct)
    assert next(rows,None) is None
    return witnesses


def verify(root):
    common = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))
    up = run_path(str(root/"ch39-rank-two-generator-up-quark.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    spec = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    weighted = run_path(str(root/"ch39-common-unitary-quark-right-trace.py"))["weighted"]
    n = 2_306_888_100_086
    _, label, _ = common["verify_nodes"](root,((n,".0294"),),".1",(),((6,7,10,13),(6,7,10,12)))
    ctx.prec = 640
    source, target = [spec(shell,-4,k)[1] for k in (n,n+1)]
    h = weighted(source,label[n,-4])
    four = acb_mat([[target[i,j] for j in range(19)] for i in (14,8,7,6)])
    eye = up["eye"](19)
    center = arb(".00000615")
    for group, transport_adj in zip(BLOCKS,build_inputs(root),strict=True):
        assert (up["adj"](transport_adj)*transport_adj).contains(eye)
        transported = transport_adj*target
        energy = up["adj"](transported)*h*transported
        amplitudes = [energy[j,j].real for j in range(19)]
        assert min(amplitudes) > 0
        raw = weighted(four,amplitudes)
        values = sorted((v.real for v in raw.eig(algorithm="rump")),key=lambda v:v.mid())
        assert all(a < b for a,b in zip(values,values[1:]))
        assert values[0] < arb(".00000601") < center < arb(".00000629") < values[1]
        print("保持的站位与四站算子升序本征值",group,[v.str(15) for v in values],flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
