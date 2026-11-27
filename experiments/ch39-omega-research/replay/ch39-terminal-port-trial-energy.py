"""五种荷各分两块覆盖双钟环面，验证三维试探空间的能量与八站读数。"""
from concurrent.futures import ProcessPoolExecutor
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ROOT = Path(__file__).resolve().parent
QUARTER = 24
CHARGES = (1, -4, 2, -3, 6)


def restricted_lower(m, sine):
    return m*(1-sine**2).sqrt() - (1-m**2).sqrt()*sine


def lower_millionths(value):
    return int((value.lower()*10**6).floor().unique_fmpz())


def circle_cells():
    def cayley(t):
        return acb(1-t*t, 2*t)/(1+t*t)

    cells = []
    for quadrant in range(4):
        for k in range(QUARTER):
            midpoint = cayley(arb(2*k+1)/(2*QUARTER))
            ends = [cayley(arb(n)/QUARTER) for n in (k, k+1)]
            # 一个象限内，格心到点的弦长在两侧各随圆弧距离递增。
            radius = max(abs(midpoint-z).upper() for z in ends)
            cells.append((acb(0,1)**quadrant * midpoint, radius))
    return cells


def certify_half(case):
    charge, half, sites = case
    graph = run_path(str(ROOT / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(ROOT / "ch39-native-light-projector-pump.py"))
    least_singular = run_path(str(ROOT / "ch39-terminal-port-light-bands.py"))["least_singular"]
    ctx.prec = 144
    adjoint = pump["adjoint"]
    triple = [sites.index(i) for i in (0, 3, 5)]
    base, links = pump["clock_matrices"](graph, charge)
    cells = circle_cells()
    sectors = []
    for frame in pump["rotation_bases"](graph):
        restricted = acb_mat([[frame[i,k] for k in range(frame.ncols())] for i in sites])
        blocks = [adjoint(frame)*matrix*frame for matrix in links]
        clocks = [[block*z + adjoint(block*z) for z, _ in cells] for block in blocks]
        sectors.append((restricted, adjoint(frame)*base*frame, *clocks))

    minima = None
    for i in range(half*len(cells)//2, (half+1)*len(cells)//2):
        for j in range(len(cells)):
            change = arb(2).sqrt()*(cells[i][1]+cells[j][1]) + 6*arb("1e-8")
            columns, gaps, light, heavy = [], [], [], []
            for restricted, fixed, first, second in sectors:
                matrix = fixed + first[i] + second[j]
                vector, gap = pump["first_line"](matrix)
                column = acb_mat([[z] for z in vector])
                energy = (adjoint(column)*matrix*column)[0,0].real
                columns.append((restricted*column).entries())
                gaps.append(gap)
                light.append(energy.upper())
                heavy.append((energy+gap).lower())
            # 所有背景在三个字符中各有一条轻线。Weyl 界覆盖整个格以及软辐边。
            band_gap = min(heavy) - max(light) - 2*change
            assert min(gaps) > 2*change and band_gap > 0
            rows = acb_mat(columns).transpose()
            center = least_singular(rows)
            tilt = change/(min(gaps)-change)
            m = restricted_lower(center, tilt).lower()
            assert 0 < m < 1 and tilt < 1
            site_rows = rows.tolist()
            triple_center = least_singular(acb_mat([site_rows[k] for k in triple]))
            triple_bound = restricted_lower(triple_center, tilt)
            # 秩三的总重带权重至多为平均能量超额的三倍除以带隙。
            critical = band_gap*m**2/3
            at_one_cent = restricted_lower(m, (arb(".03")/band_gap).sqrt())
            bounds = tuple(map(lower_millionths, (m, band_gap, critical, at_one_cent, triple_bound)))
            minima = bounds if minima is None else tuple(map(min, zip(minima, bounds, strict=True)))
    return charge, half, minima


def verify():
    sites = run_path(str(ROOT / "ch39-terminal-port-support.py"))["verify"](ROOT)
    cases = [(q, h, sites) for q, h in product(CHARGES, range(2))]
    with ProcessPoolExecutor(max_workers=len(cases)) as pool:
        results = list(pool.map(certify_half, cases))
    for charge, half, values in results:
        print("荷、半环编号、八站奇异值/轻重隙/失去单射所需超额/超额0.01时八站奇异值/三站奇异值（百万分之一）",
              charge, half, *values)
    joint = tuple(map(min, zip(*(row[2] for row in results), strict=True)))
    assert joint[0] > 478280 and joint[1] > 1387850 and joint[2] > 112690 and joint[3] > 348970
    assert joint[4] > 297470
    print("全环面共同下界（百万分之一）", *joint)


if __name__ == "__main__":
    verify()
