"""一个整数线性泛函将 1,136 个子壳类顶点与九个十二次目标的整个凸包严格分开。"""
from itertools import product
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat, ctx


# 仅使用 x=t/9-1 的一至十一次幂。首项统一归一化，常数项自由。
NORMAL = (-4294967296, 1754467993, -411676551, -145440381, 218247272,
          -81477457, -49221441, 88163752, -47067236, -49971991, 189104479)


def verify(replay):
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))

    sector = load("geometric-selector-mixture")
    grouped = load("geometric-exchangeable-priors")
    shell = load("three-ray-chiral-boundary")
    pump = load("native-light-projector-pump")
    light = load("light-spectrum-zero-selector")["light_spectra"](replay)
    ctx.prec = 320
    triple = sector["geometry"](shell,pump)
    factors = sector["squared_factors"]
    labels = list(product(grouped["MASKS"],grouped["MODES"]))
    row = {pair:i for i,pair in enumerate(labels)}
    charges = []
    for charge in (-4,2):
        polynomials = [factors(triple(charge,mask,mode)) for mask,mode in labels]
        charges.append(arb_mat([[poly[k] for k in range(7)] for poly in polynomials]))
    kernel = arb_mat([[NORMAL[i+j-1] if 1 <= i+j <= 11 else 0
                       for j in range(7)] for i in range(7)])
    individual = charges[0]*kernel*charges[1].transpose()
    vertices = []
    for pairs in grouped["geometry_classes"]().values():
        for mode, other in grouped["MODE_PAIRS"]:
            vertices.append(sum(individual[row[u,mode],row[d,other]] for u,d in pairs)/len(pairs))
    assert len(vertices) == 1136
    high = max(z.upper() for z in vertices)

    den = 10**46
    roots = {}
    for q,number in ((-4,2345666),(2,964890)):
        for shift in (-1,0,1):
            roots[q,shift] = []
            for level in light(q,number+shift):
                numerator = (den*level+arb(1)/2).floor().unique_fmpz()
                assert numerator is not None
                roots[q,shift].append(arb(numerator)/den)
    targets = []
    for i,j in product((-1,0,1),repeat=2):
        poly = factors(roots[-4,i]+roots[2,j])
        targets.append(sum(coeff*poly[k] for k,coeff in enumerate(NORMAL,1)))
    low = min(z.lower() for z in targets)
    assert high < 363504000 and low > 363546200
    assert low-high > 42200
    print("1,136 个顶点的泛函上界、九个目标的泛函下界及最小间隔",
          high.str(15),low.str(15),(low-high).str(15))


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
