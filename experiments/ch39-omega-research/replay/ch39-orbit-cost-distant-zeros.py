"""留出窗外实际反例的最优费用，及含任意高零点对的开相位邻域。"""
import json
from pathlib import Path
from runpy import run_path

from flint import acb, arb


def verify():
    replay = Path(__file__).resolve().parent
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    experiment = run_path(str(replay/"ch39-minimum-orbit-specific-probability.py"))
    optimal = experiment["setup"]()
    data = json.loads((replay/"ch39-nine-orbit-cost-bases.json").read_text())
    basis = next(row["support"] for row in data if row["offsets"] == [0,0])
    central = optimal((0,0),basis)
    offsets = (0,8244-964890)
    far = optimal(offsets,basis)
    assert arb("0.00016501496") < far < arb("0.00016501497") < central

    radius = arb("1e-11")
    phase_radius = arb("1e-13")
    robust = optimal(offsets,basis,radius)
    assert robust < arb("0.000165016") < central
    offset = shell["soft"]["point"]["study"]["offset"]
    kappa = arb(23)/25*arb(2)**-20
    # 每条钟边的单位相位变动至多 phase_radius，任一站的度不超过 18。
    # 三条软辐边所在的星的拉普拉斯绝对行和不超过 6。
    # 取两边实际高度均大于 10,000 的新增配对，允许两部门各自逼近参考角。
    for n in (2345666,8244):
        delta0 = kappa/(acb.zeta_zero(n).imag**2+offset)
        delta_max = kappa/(arb(10000)**2+offset)
        assert 0 < delta0 and 0 < delta_max
        assert 18*phase_radius+6*(delta0+delta_max)+arb("1e-46") < radius
    print("中央、远处配对的最小费用与开放邻域包围", central.str(16), far.str(16), robust.str(10))


if __name__ == "__main__":
    verify()
