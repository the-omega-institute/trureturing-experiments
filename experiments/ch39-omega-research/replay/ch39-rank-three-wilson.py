"""一条固定原边的电导改动在三个轻态上诱导不交换的两条闭路。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def frame(matrix):
    values, vectors = matrix.eig(right=True,algorithm="rump")
    ix = sorted(range(len(values)),key=lambda i:float(values[i].real))
    for i,j in zip(ix[:3],ix[1:4]):
        assert values[i].real < values[j].real
    cols = []
    for j in ix[:3]:
        n = sum(abs(vectors[i,j])**2 for i in range(19)).sqrt()
        cols.append([vectors[i,j]/n for i in range(19)])
    return acb_mat(cols).transpose()


def polar(matrix):
    gram = adj(matrix)*matrix
    values, vectors = gram.eig(right=True,algorithm="rump")
    assert all(z.imag.contains(0) and z.real>0 for z in values)
    inverse_root = acb_mat(3,3)
    for k,v in enumerate(values):
        col = acb_mat([[vectors[j,k]] for j in range(3)])
        inverse_root += (col*adj(col))/((adj(col)*col)[0,0].real*v.real.sqrt())
        assert all(not v.overlaps(z) for z in values[:k])
    return matrix*inverse_root


def measure(shell, retune):
    frames=[]
    for n in (22,114,70,18):
        mat = shell["laplacian"](-3, acb.zeta_zero(n).imag)
        assert abs(mat[0,2]).contains(1)
        mat[0,0] += retune
        mat[2,2] += retune
        mat[0,2] *= 1+retune
        mat[2,0] = mat[0,2].conjugate()
        frames.append(frame(mat))
    def edge(i,j):
        return polar(adj(frames[j])*frames[i])
    w1 = edge(2,0)*edge(1,2)*edge(0,1)
    w2 = edge(3,0)*edge(1,3)*edge(0,1)
    diff = w1*w2-w2*w1
    return sum(abs(z)**2 for z in diff.entries()).nonnegative_part().sqrt()


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 224
    without = measure(shell, arb(0))
    with_edge = measure(shell, arb(1)/8)
    assert without.is_finite() and without.contains(0) and without < arb("1e-30")
    assert arb("5.95590e-6") < with_edge < arb("5.95592e-6")
    print("荷 −3，闭路 (22,114,70,22) 与 (22,114,18,22)")
    print("原壳及单边电导增加 1/8 后的 Wilson 幺正对易子距恒等 Frobenius 范数",
          without.str(12),with_edge.str(18))


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
