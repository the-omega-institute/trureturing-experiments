"""在 ±64 邻域所有能重现六质量的配对中，比较冻结的轨道身份费用 η。

用原六质量分类重新检验全部 16,641 对。对保留的 97 对，分别以各自的
六个实际能量取整后校准目标多项式。每对均以原始可行混合和可行对偶
包围最优费用，两种证明使用同一组十二列。排序只使用严格分离的费用区间。
"""
from pathlib import Path
from runpy import run_path
import json


REPLAY = Path(__file__).resolve().parent


def basis_map(feasible):
    bases = {}
    # 九个近邻与十六个外圈已有公开基，其中两对通过六质量筛选。
    for filename in ("ch39-nine-orbit-cost-bases.json", "ch39-holdout-orbit-cost-bases.json"):
        for record in json.loads((REPLAY / filename).read_text()):
            offsets = tuple(record["offsets"])
            if offsets in feasible:
                bases[offsets] = record["support"]
    assert set(bases) == {(0, 0), (0, -2)}

    # 若多个目标位于同一原始单形，数据仅保存一次列清单。
    groups = json.loads((REPLAY / "ch39-mass-feasible-orbit-cost-bases.json").read_text())
    for group in groups:
        for pair in group["offsets"]:
            offsets = tuple(pair)
            assert offsets not in bases
            bases[offsets] = group["support"]
    assert set(bases) == feasible
    return bases


def verify():
    classify = run_path(str(REPLAY / "ch39-zero-pair-one-mixing-angle.py"))["verify"]
    exact, added = classify()
    feasible = set(exact) | set(added)
    assert len(exact) == 93 and len(added) == 4 and len(feasible) == 97
    assert all(up == 0 for up, down in feasible)
    bases = basis_map(feasible)
    optimal = run_path(str(REPLAY / "ch39-minimum-orbit-specific-probability.py"))["setup"]()
    costs = {offsets: optimal(offsets, support) for offsets, support in bases.items()}
    # 中点只提出顺序；相邻区间严格分离才确定全部名次。
    ranked = sorted(costs, key=lambda k: float(costs[k].mid()))
    assert all(costs[a] < costs[b] for a, b in zip(ranked, ranked[1:]))
    assert ranked[:3] == [(0, -58), (0, 0), (0, 39)]
    assert ranked[0] in exact and ranked[1] in exact
    for rank, pair in enumerate(ranked, 1):
        status = "六质量中心" if pair in exact else "仅完整质量窗"
        print("名次、编号偏移、质量条件、最小具名份额", rank, pair, status, costs[pair].str(21))


if __name__ == "__main__":
    verify()
