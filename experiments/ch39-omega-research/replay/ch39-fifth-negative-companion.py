"""第五壳正动能与弱差要求负权层及三条强负权层分别进入低能区。"""
from pathlib import Path
from runpy import run_path

from flint import arb
import numpy as np


def verify(root):
    core = run_path(str(root / "ch39-fourth-third-pivot-global-kinetic-obstruction.py"))
    energies = (1778, 3162, 10000, 31623, 100000)
    lower, upper, _, weights, initial, _ = core["build_grid"](
        root, (0,), number=2_306_888_100_090, energies=energies)
    ticks = core["cutoff_ticks"]()
    scale = core["SCALE"]
    whole = (0, len(ticks)-1)
    neg = [i for i,w in enumerate(weights) if w < 0]
    assert [i+1 for i,w in enumerate(weights) if w > 0] == [1, 4, 10, 16]
    # 用精确阈值左右两侧的网格点包住最后一段弦及其舍入余量。
    def split(gev):
        tick = 1024*(arb(gev)/arb("91.1876")).log()
        if tick < int(ticks[0]):
            return 0,0
        if tick > int(ticks[-1]):
            return whole[1],whole[1]
        index = int(np.searchsorted(ticks, float(tick.mid())))
        assert int(ticks[index-1]) < tick < int(ticks[index])
        return index-1, index

    def band_upper(column, bands):
        result = initial[column]
        for i,(left,right) in enumerate(bands):
            values = slice(left, right+1)
            tip = (int(upper[values,column].max()) if weights[i]>0
                   else int(lower[values,column].min()))
            result += weights[i]*tip/scale
        return result.upper()

    for cap, locked, positives in [
        (4000, neg, ((0,1150,3162),(9,350,1778),(15,2800,10000))),
        (30000, [16,17,18], ((0,4000,31623),(9,260,10000),(15,12000,100000))),
    ]:
        bands = [whole for _ in weights]
        for i in locked:
            bands[i] = (split(cap)[0],whole[1])
        for i, threshold, energy in positives:
            varied = bands.copy()
            varied[i] = (0,split(threshold)[1])
            column = 4+2*energies.index(energy)
            bound = 10*band_upper(column,varied)
            print("负层假设、正层、其截止能量上限、测试能标、色读数上界",
                  cap,i+1,threshold,energy,bound.str(11),flush=True)
            assert bound < 0
        # 持续正色迫使三个正权层分别高于各自阈值。
        for i,threshold,_ in positives:
            bands[i] = (split(threshold)[0],whole[1])
        # 五式与持续正动能的既有引力必要条件给出 x_4>0.10，此处仅取更弱的下界。
        bands[3] = (split(10**18)[0],whole[1])
        weak = band_upper(2,bands)/100
        print("三层均高时弱差上界", cap, weak.str(11), flush=True)
        assert weak < 0

    for i, threshold, locked, cap, energy in [
        (0, 493, [16,17,18], 2800, 1778),
        (9, 249, [16,17,18], 33000, 10000),
        (15, 1025, [17], 3400, 10000),
    ]:
        bands = [whole for _ in weights]
        bands[i] = (0,split(threshold)[1])
        for j in locked:
            bands[j] = (split(cap)[0],whole[1])
        bound = 10*band_upper(4+2*energies.index(energy), bands)
        print("指定正层和截止、指定负层和截止、测试能标、色读数上界",
              i+1,threshold,[j+1 for j in locked],cap,energy,bound.str(12),flush=True)
        assert bound < 0


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
