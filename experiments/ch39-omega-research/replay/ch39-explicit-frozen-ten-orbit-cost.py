"""实际零点对的冻结十项读出与具名轨道最小费用。

质量、三只 CKM 角及 CP 相位使用同一复源 K、同一制备算子及原部门标尺。
费用沿用各零点对自身六个能量取整至小数后 46 位的首一校准多项式。
"""
from functools import cache
import json
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


CASES = (
    # 上、下型零点编号，逐项最大偏差/参考半宽的上界，费用改进下界。
    (27804036, 964890, ".101", "9.46e-8"),
    (27804036, 57113204, ".199", "1.51e-7"),
)


def verify():
    replay = Path(__file__).resolve().parent
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))

    model = load("two-zero-quark-window")
    at_zero = load("quadratic-gram-source")["at_zero"]
    optimal = load("minimum-orbit-specific-probability")["setup"]()
    candidates = json.loads((replay / "ch39-nine-orbit-cost-bases.json").read_text())
    basis = next(v["support"] for v in candidates if v["offsets"] == [0, 0])
    center = optimal((0, 0), basis)
    scales = {q: arb(scale) for q, _, _, scale in model["PHASE_BOXES"]}
    ctx.prec = 384

    @cache
    def spectrum(charge, index):
        # 在该实际零点取三个最低字符本征投影，沿用原源列与谱滤波。
        prepared = at_zero(model, charge, index)
        yukawa = prepared * model["K"] * scales[charge]
        return model["eig"](model["adj"](yukawa) * yukawa)

    for up, down, error_bound, gain_bound in CASES:
        u, left = spectrum(-4, up)
        d, right = spectrum(2, down)
        values = [v.sqrt() for v in u+d] + model["mixing_angles"](left, right)
        deviations = [abs(v-target)/width for v, target, width
                      in zip(values, model["ref"], model["tol"], strict=True)]
        assert len(deviations) == 10 and all(error < arb(error_bound) for error in deviations)

        # 原始权重为正，冻结的十二列同时给出可行对偶下界。
        cost = optimal((up-2345666, down-964890), basis)
        gain = center-cost
        assert gain > arb(gain_bound)
        print("上、下型编号与十项归一偏差", up, down, *(x.str(10) for x in deviations))
        print("最大偏差上界、最小具名份额及相对中央的费用改进",
              max(x.upper() for x in deviations).str(14), cost.str(18), gain.str(14))


if __name__ == "__main__":
    verify()
