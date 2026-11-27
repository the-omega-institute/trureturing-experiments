"""仅数值：64x64 联合相位格点上的跨荷轻空间读数（前三模定义，delta = 0，浮点）。

这些读数不是证书，只用于说明中标为“仅数值”的几句话。
"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parent
SHELL = run_path(str(ROOT / "ch39-three-ray-chiral-boundary.py"))
GEO, SOFT = SHELL["geo"], SHELL["soft"]
SITES, POLE = GEO["sites"], SOFT["pole"]
EDGES = [(i, j, d, complex(z.mid())) for i, j, d, z in GEO["edges"]]
CHARGES = (1, 2, -3, -4, 6)
ODD_EVEN = [(1, 2), (1, -4), (1, 6), (2, -3), (-3, -4), (-3, 6)]
N = 64
ROT = [SITES.index((y, z, x)) for x, y, z in SITES]
INV = [ROT.index(j) for j in range(19)]                 # (R x)[j] = x[R^{-1} j]


def lap(q, X, Y):
    s = abs(q)
    M = np.zeros((19, 19), complex)
    for i, j, d, z in EDGES:
        M[i, i] += 1
        M[j, j] += 1
        ph = 1 if POLE in (i, j) else z**s*(np.exp(1j*s*X) if d == 724 else np.exp(1j*s*Y) if d == 804 else 1)
        M[i, j] -= ph
        M[j, i] -= np.conj(ph)
    return M if q > 0 else M.conj()


def light(q, X, Y):
    w, v = np.linalg.eigh(lap(q, X, Y))
    o = np.argsort(np.abs(w*(w-15)))[:3]
    return v[:, o]


def sector(x):
    r = np.vdot(x, x[INV])
    return int(round(np.angle(r)/(2*np.pi/3))) % 3


def main():
    leak = {p: np.zeros((N, N)) for p in combinations(CHARGES, 2)}
    phase, mono = {p: np.zeros((N, N)) for p in ODD_EVEN}, {p: np.zeros((N, N)) for p in ODD_EVEN}
    sec0 = 0.0
    for a in range(N):
        for b in range(N):
            X, Y = 2*np.pi*a/N, 2*np.pi*b/N
            V = {q: light(q, X, Y) for q in CHARGES}
            P = {q: V[q]@V[q].conj().T for q in CHARGES}
            for p in leak:
                q, r = p
                leak[p][a, b] = 1-np.trace(P[q]@P[r]).real/3
            for q, r in ODD_EVEN:
                A, B = np.abs(P[q]), np.abs(P[r])
                phase[q, r][a, b] = 1-np.sum(A*B.T)/3
                S = -np.sort(-A, axis=1)@(-np.sort(-B, axis=1)).T
                i, j = linear_sum_assignment(-S)
                mono[q, r][a, b] = 1-S[i, j].sum()/3
                x0 = [V[q][:, k] for k in range(3) if sector(V[q][:, k]) == 0]
                y0 = [V[r][:, k] for k in range(3) if sector(V[r][:, k]) == 0]
                assert len(x0) == len(y0) == 1
                sec0 = max(sec0, abs(np.vdot(x0[0], y0[0]))**2)
    pmin = min(ODD_EVEN, key=lambda p: phase[p].min())
    mmin = min(ODD_EVEN, key=lambda p: mono[p].min())
    print("调相下界泛函在奇偶六对上的格点最小值", round(phase[pmin].min(), 4), pmin)
    print("调相—重排下界泛函在奇偶六对上的格点最小值", round(mono[mmin].min(), 4), mmin,
          np.unravel_index(mono[mmin].argmin(), (N, N)))
    print("奇偶六对扇区 0 的 cos^2 格点最大值", round(sec0, 4))
    argmins = {p: tuple(int(v) for v in np.unravel_index(leak[p].argmin(), (N, N))) for p in ODD_EVEN}
    print("奇偶六对泄漏的格点最小点", argmins)
    assert 0.30 < phase[pmin].min() < 0.32 and pmin == (-3, 6)
    assert 0.19 < mono[mmin].min() < 0.21 and mmin == (-3, 6)
    assert sec0 < 0.03
    assert all(argmins[p] == (N//2, N//2) for p in ODD_EVEN if p != (-3, 6)) and argmins[-3, 6] != (N//2, N//2)


if __name__ == "__main__":
    main()
