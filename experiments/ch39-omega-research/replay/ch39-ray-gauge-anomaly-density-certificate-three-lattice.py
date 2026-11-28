"""3^4 周期格上纯超荷凸片 K_a 内一点的十五分量手征密度 A 的 Arb 区间证书。

复用 ch39-ray-gauge-anomaly-density-certificate.py 的全部函数（同一 A 的定义、同一 Bauer-Fike 谱区间与
Legendre 符号函数展开）。场 A = T/5000，T 为下方逐字给出的 +-1 模式，
证书格点为 30 号格点（字典序坐标 (1,0,1,0)）。精度 700 位，截断误差 <= 1e-25。
"""
import os
import runpy
import time
from fractions import Fraction as Fr

from flint import arb, ctx

M = runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), "ch39-ray-gauge-anomaly-density-certificate.py"),
                   run_name="ray_gauge_anomaly_module")
LS = (3, 3, 3, 3)
SITE = 30
PATTERN = ("+++++++----++--+-+--+---+-+--+++++++++-+--+-++-----++-++++-+----+---+++-+---+++-+-++-+++++-+---++++++-++--+-+----+--++-+"
           "-+-++++-+-----+--++--+---+-++-++-+++++--+-++-+---+++++--+---+-+++-++---+--+-+-+--+-++--+-++--+-++-+++++++--+--+++++-+-++"
           "+-++-+++-+----+------++++----+---+-++++++++++-++--++-----++++-++++-+-+++--+-+------+")


def run():
    t0 = time.time()
    ctx.prec = 700
    M["TOL_LOG10"] = -25
    M["sign_diagonal"].__globals__["TOL_LOG10"] = -25
    V, nb = M["lattice"](LS)
    assert V == 81 and len(PATTERN) == 4*V
    field = M["field_from_pattern"](LS, PATTERN, M["FIELD_SCALE"])
    assert max(abs(v) for v in field.values()) == Fr(1, 5000) < Fr(1, 4800)
    pm = M["plaquette_max"](LS, field)
    assert 6*pm <= Fr(1, 200)
    log = []
    A, gap = M["anomaly_density"](LS, field, [SITE], log)
    A0 = A[0]
    Nmax = max(e[6] for e in log)
    print(f"[格 {LS}，V={V}] 格点 {SITE}；最大方格角 {pm}；min|mu| >= {gap.lower().str(8)}；Legendre 项数 <= {Nmax}")
    print(f"  A({SITE}) = {A0.str(15)}")
    assert pm == Fr(1, 1250) and gap > arb("0.999086") and Nmax <= 155
    assert gap > arb(199)/200
    assert A0.rad() < arb(10)**-20
    assert A0 < arb("-9.0757234e-13") and A0 > arb("-9.0757235e-13")
    print(f"耗时 {time.time()-t0:.1f} s")


if __name__ == "__main__":
    run()
    print("全部断言通过")
