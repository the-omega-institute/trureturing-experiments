"""固定同一组十二种原壳配方，仅调整概率即可唯一选择九种邻号配对。"""
from functools import cache
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


SUPPORT = (
    (84, 72, "D"), (2, 2, "I"), (42, 52, "I"), (31, 92, "I"),
    (1, 60, "I"), (60, 60, "I"), (60, 123, "I"), (62, 124, "I"),
    (63, 127, "D"), (92, 90, "I"), (111, 127, "D"), (120, 2, "I"),
)


def verify():
    replay = Path(__file__).resolve().parent
    mixture = run_path(str(replay/"ch39-geometric-selector-mixture.py"))
    shell = run_path(str(replay/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay/"ch39-native-light-projector-pump.py"))
    read = run_path(str(replay/"ch39-light-spectrum-zero-selector.py"))
    light = read["light_spectra"](replay)
    # 与九对唯一性证书共享十进制取整的有理数；随后提升运算精度。
    @cache
    def rational_roots(charge, index):
        numerators = [(v*10**46+arb(1)/2).floor().unique_fmpz() for v in light(charge,index)]
        assert all(n is not None for n in numerators)
        return numerators

    centers = ((-4,2345666), (2,964890))
    for charge, center in centers:
        for offset in (-1,0,1):
            rational_roots(charge, center+offset)
    ctx.prec = 640
    triple = mixture["geometry"](shell, pump)
    squared_factors = mixture["squared_factors"]
    recipes = [squared_factors(triple(-4,u,m)+triple(2,d,m)) for u,d,m in SUPPORT]
    inverse = mixture["simplex_inverse"](recipes)
    for du, dd in product((-1,0,1), repeat=2):
        roots = rational_roots(-4,2345666+du)+rational_roots(2,964890+dd)
        target = squared_factors([arb(n)/10**46 for n in roots])
        weights, smallest, radius = mixture["barycentric"](inverse, target)
        _, error = mixture["rationalize"](weights, recipes, target)
        # 九个选择器对全体实际编号的作用差均大于 8.11e-71。
        assert smallest > arb("2.37e-8") and error < arb("8.11e-71")
        print("编号偏移、最小权、到单形面的距离、权重有理化误差",
              (du,dd), smallest.str(9), radius.str(9), error.str(9))


if __name__ == "__main__":
    verify()
