"""检验同一六根取整规则能分别唯一选中上下各三个相邻编号的任意配对。"""
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    replay = Path(__file__).resolve().parent
    selector = run_path(str(replay/"ch39-rational-unique-light-action.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    read = run_path(str(replay/"ch39-light-spectrum-zero-selector.py"))
    light = read["light_spectra"](replay)
    ctx.prec = 256
    denominator = 10**46
    eps = arb("1e-32")

    @lru_cache(maxsize=None)
    def calibrated_roots(charge, number):
        levels = light(charge, number)
        roots = []
        for level in levels:
            numerator = (denominator*level + arb(1)/2).floor().unique_fmpz()
            assert numerator is not None
            roots.append(arb(numerator)/denominator)
        return levels, roots

    for du, dd in product((-1,0,1), repeat=2):
        centers = ((-4,2345666+du), (2,964890+dd))
        actual, roots = [], []
        for charge,number in centers:
            energy, calibrated = calibrated_roots(charge, number)
            actual.extend(energy)
            roots.extend(calibrated)
        target_intervals = [r+arb(0,eps.upper()) for r in roots]
        selector["verify_number_boxes"](pump, shell, centers, target_intervals)
        # 若第 i 根最近，则至其他每根的距离至少是对应根间距之半。
        lower = eps**2 * min(
            prod((abs(x-y)/2)**2 for j,y in enumerate(roots) if i!=j).lower()
            for i,x in enumerate(roots)
        )
        upper = sum(prod((e-r)**2 for r in roots) for e in actual)
        assert lower-upper > arb("8.11e-71")
        print("上下型校准编号与统一作用差下界", centers, (lower-upper).str(12))


if __name__ == "__main__":
    verify()
