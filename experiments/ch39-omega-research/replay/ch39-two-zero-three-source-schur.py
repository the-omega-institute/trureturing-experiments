"""用六质量的归一化迹约束逐对检验三源相干幅度与相位。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb


def reduce_pair(form, coordinates, adjoint, i, j, k):
    # 减去第 k 列的倍数不改变其 Schur 补，有助于分开原三源的近平行列。
    initial = adjoint(coordinates)*form*coordinates
    assert initial[k,k].real > 0
    shifted = acb_mat(coordinates)
    for col in (i,j):
        multiple = initial[k,col].mid()/initial[k,k].real.mid()
        for row in range(3):
            shifted[row,col] -= multiple*coordinates[row,k]
    matrix = adjoint(shifted)*form*shifted
    pivot = matrix[k,k].real
    assert pivot > 0
    return acb_mat([[matrix[p,q]-matrix[p,k]*matrix[k,q]/pivot
                     for q in (i,j)] for p in (i,j)])


def verify():
    inputs = run_path(str(Path(__file__).with_name("ch39-two-zero-source-correlations.py")))
    (up,down), coordinates, factor, adj = inputs["certified_inputs"]()
    form = arb("74.25")*adj(down)*down-factor*adj(up)*up
    rows = (
        ("正交源", (".93946",".97632",".98253"), (".350",".219",".188"), ".9322"),
        ("原算术源", (".99980",".99908",".99972"), (".020",".043",".024"), ".99906"),
    )
    for (name, lower_bounds, half_angles, triple_bound), transform in zip(rows,coordinates,strict=True):
        for k, ((i,j),lower,half_angle) in enumerate(zip(
                ((1,2),(2,0),(0,1)),lower_bounds,half_angles,strict=True)):
            reduced = reduce_pair(form,transform,adj,i,j,k)
            a, b, z = reduced[0,0].real, reduced[1,1].real, reduced[0,1]
            assert a > 0 and b > 0 and abs(z) > 0
            projected = (a*b).sqrt()/abs(z)
            phase = (-z).arg()
            assert arb(lower) < projected < 1
            assert projected.acos() < arb(half_angle)
            assert projected > abs(phase.cos())
            print(name,i+1,j+1,"相关在指定相位上的投影下界、相位中心与半宽",
                  projected.str(9),phase.str(10),projected.acos().str(8))
        assert 2*sum(map(arb,lower_bounds))/3-1 > arb(triple_bound)


if __name__ == "__main__":
    verify()
