"""六个实际高度的六个重场部门具有可分辨的谱线及满秩的规范权。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat


def verify():
    root = Path(__file__).resolve().parent
    transfer = run_path(str(root/"ch39-five-zero-gravity-gauge-transfer.py"))
    indices = (transfer["BASE"], *transfer["CANDIDATES"], transfer["BASE"]+1)
    spectra = run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))
    _, lib, rows = spectra["gravity_fields"](root, indices, tight_neutral=True)
    model, _, _, _, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    coefficients = run_path(str(root/"ch39-signed-gravity-unification.py"))["gauge_weights"](model)
    blocks = [(start, start+16) for start in range(0, 80, 16)] + [(80, 98)]
    groups = [row[start+1:end+1] for row in rows for start, end in blocks]
    assert len(groups) == 36

    # 每组只需一条区别于其他 35 组的正质量线。同组出现简并时，系数乘以重数。
    for index, group in enumerate(groups):
        start, end = blocks[index % len(blocks)]
        g, alpha, witness = group[0]
        assert g != 0 and not alpha.contains(0) and witness > 0
        assert all(h == g and (a-alpha).contains(0) for h, a, _ in group)
        gauge = coefficients[start]
        assert any(abs(v) > 0 for v in gauge)
        for row in coefficients[start:end]:
            assert all((arb(a)-arb(b)).contains(0)
                       for a, b in zip(row, gauge, strict=True))
        other_lines = [z for k, other in enumerate(groups) if k != index for _, _, z in other]
        assert all(not witness.overlaps(z) for z in other_lines)

    # 前三个表示的三条规范向量张成全部三支阈值。
    gauge_basis = arb_mat([coefficients[start] for start, _ in blocks[:3]])
    assert not gauge_basis.det().contains(0)
    equations = 6+3*len(blocks)*(len(indices)-1)
    assert equations == 96 and len(indices)**len(blocks) == 46_656
    print("六个编号", *indices)
    print("具有独占正谱线的部门-高度组数", len(groups), "三规范向量行列式", gauge_basis.det().str(14))
    print("共同约束的独立行数", equations, "覆盖完整重谱数", len(indices)**len(blocks))


if __name__ == "__main__":
    verify()
