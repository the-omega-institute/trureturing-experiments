"""两种四站味空间与一条全轴有界函数覆盖三枚实际零点的全部味和引力读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify(root, *, return_data=False):
    pair = run_path(str(root / "ch39-two-zero-simultaneous-quark-gravity.py"))
    right_frame, training = pair["verify"](root)
    left = 2_306_888_100_084
    flavor = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))
    geometry = (".0182", ".113", ("-.21668", "-.53516", ".71837", "-.38807"))
    left_frame, responses, _ = flavor["verify_nodes"](
        root, ((left, "0"),), ".74", geometry, ((1,10,16,18), (6,7,10,18))
    )
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ordered = run_path(str(root / "ch39-bounded-common-rational-amplitude.py"))["ordered"]
    gravity = run_path(str(root / "ch39-left-cutoff-two-turns.py"))
    levels, positions = gravity["verify"](root)
    records = training + list(zip(levels, [x.sqrt() for x in positions], strict=True))
    ctx.prec = 320
    height = acb.zeta_zero(left).imag
    for q in (-4, 2):
        records.extend(zip(ordered(shell, q, height), responses[left,q], strict=True))
    records.sort(key=lambda row: row[0].mid())
    floor = arb("91.1876") / arb("1.22089e19")
    assert len(records) == 171
    assert all(y-x > arb("8.28787e-7") for (x,_),(y,_) in zip(records,records[1:]))
    assert all(floor < v < arb(".999") for _,v in records)
    print("三个编号的校准节点数及最大振幅", len(records),
          max(v.upper() for _,v in records).str(14), flush=True)
    if return_data:
        return records, (left_frame, right_frame)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
