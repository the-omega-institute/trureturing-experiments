"""任意十四站上的两种层级，以及任意十八站对上下夸克矩阵的满射。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat,arb,ctx


def verify():
    root = Path(__file__).parent
    flavor = run_path(str(root/"ch39-minimum-flavor-support.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 256
    up,down = flavor["flavor_maps"](shell,11431)
    matrix = acb_mat(up.tolist()+down.tolist())
    minors = []
    for remove in range(19):
        square = acb_mat([[matrix[i,j] for j in range(19) if j != remove] for i in range(18)])
        minors.append(abs(square.det()))
    assert min(minors) > arb("1.9499e-31")
    assert max(minors) < arb("1.7556e-29")
    values = (matrix*flavor["adjoint"](matrix)).eig()
    least = min(v.real for v in values).sqrt()
    greatest = max(v.real for v in values).sqrt()
    assert arb(".00011407") < least < arb(".00011408")
    assert arb(".73731") < greatest < arb(".73732")
    print("上下型联合矩阵最小十八阶余子式、最小和最大奇异值",
          min(minors).str(12),least.str(12),greatest.str(12))


if __name__ == "__main__":
    verify()
