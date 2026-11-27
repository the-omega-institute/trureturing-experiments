"""同一十二站支撑逐层实现十三站见证的十九个目的振幅，逐层容许不同酉矩阵。"""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb


SITES = (0, 1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 18)


def direction(block, coupling, radius):
    # 浮点球面二次型计算只选方向；最终严格不等式用归一化后的高斯整数向量核验。
    eigen, basis = np.linalg.eigh(block)
    b = basis.conj().T @ coupling
    lo = -eigen[0]
    hi = np.linalg.norm(b) / radius + 1e-8
    for _ in range(56):
        mid = (lo + hi) / 2
        if np.linalg.norm(b / np.maximum(eigen + mid, 1e-21)) > radius:
            lo = mid
        else:
            hi = mid
    x = -b / np.maximum(eigen + hi, 1e-21)
    missing = max(0, radius**2 - np.vdot(x, x).real)
    phase = 1j*x[0]/abs(x[0]) if abs(x[0]) > 1e-21 else 1
    x[0] += missing**.5 * phase
    vec = basis @ x
    return vec / np.linalg.norm(vec)


def verify(root):
    character = run_path(str(root / "ch39-up-unrestricted-character-floor.py"))
    original = run_path(str(root / "ch39-local-unitary-up-mass-witness.py"))
    up = run_path(str(root / "ch39-rank-two-generator-up-quark.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spectral = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    weights, target = character["shell_data"](root)
    _, source = spectral["spectrum"](shell, -4, up["N"] + 1)
    adj = up["adj"]
    h = acb_mat([[source[i, j] * weights[j] for j in range(19)] for i in range(19)]) * adj(source)
    transport = original["local_transport"](adj, up["eye"])
    goals = up["read"](source, target, weights, transport)
    h_mid = np.array(h.tolist(), dtype=complex)
    block = h_mid[np.ix_(SITES, SITES)]
    below, above = [], []
    for j, goal in enumerate(goals):
        native = acb_mat([[target[i, j]] for i in range(19)])
        outside = acb_mat(native)
        for i in SITES:
            outside[i, 0] = 0
        radius = sum(abs(native[i, 0])**2 for i in SITES).sqrt()
        assert radius > 0
        b = (h_mid @ np.array(outside.entries(), dtype=complex))[list(SITES)]
        nominal = (adj(native) * h * native)[0, 0].real
        for sign, margins in ((1, below), (-1, above)):
            if sign*(goal - nominal) > 0:
                value = nominal
            else:
                vector = direction(sign*block, sign*b, float(radius.mid()))
                seed = acb_mat([[acb(int(round(z.real*1e10)), int(round(z.imag*1e10)))] for z in vector])
                inner = seed * (radius / (adj(seed) * seed)[0, 0].real.sqrt())
                full = acb_mat(outside)
                for row, site in enumerate(SITES):
                    full[site, 0] = inner[row, 0]
                assert (adj(full)*full)[0, 0].contains(1)
                value = (adj(full) * h * full)[0, 0].real
            assert sign*(goal - value) > 0, (j+1, sign)
            margins.append((sign*(goal - value)).lower())
    variation = arb("7e-10")
    assert min(below) > variation and min(above) > variation
    print("活动站位", SITES)
    print("十九层上下试点到各自目标的最小严格余量", min(below).str(12), min(above).str(12))
    four = acb_mat([[target[i, j] for j in range(19)] for i in (14, 8, 7, 6)])
    raw = acb_mat([[four[i, j] * goals[j] for j in range(19)] for i in range(4)]) * adj(four)
    levels = sorted((z.real for z in raw.eig(algorithm="rump")), key=lambda x: x.mid())
    for j, (lo, hi) in enumerate(up["WINDOWS"]):
        assert levels[j] + variation < arb(lo) < arb(hi) < levels[j+1] - variation
    up["compression"](raw, list(map(arb, up["CENTERS"])))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
