"""直接采用原始 Gram 或只在每个成分中混合两条源时的夸克迹比。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    root = Path(__file__).resolve().parent
    (au,ad), coords, c, adj = run_path(str(root/"ch39-two-zero-source-correlations.py"))["certified_inputs"]()
    for change,limits in zip(coords,(("3.55","3.50","3.44"),("3.58","3.58","3.70")),strict=True):
        u,d = au*change,ad*change
        gu,gd = adj(u)*u,adj(d)*d
        for (i,j),limit in zip(((0,1),(0,2),(1,2)),limits,strict=True):
            pencil = arb(limit)*gd-c*gu
            assert pencil[i,i].real > 0
            assert pencil[i,i].real*pencil[j,j].real > abs(pencil[i,j])**2
    # 原三源 B 位于正交源张成的空间内，故 B*B = change* change（星号指伴随）。
    raw_gram = adj(change)*change
    value = c*(gu*raw_gram).trace().real/(gd*raw_gram).trace().real
    assert arb("1.0805159") < value < arb("1.0805163")
    print("原 Gram 对应的六质量迹比",value.str(10))


if __name__ == "__main__":
    verify()
