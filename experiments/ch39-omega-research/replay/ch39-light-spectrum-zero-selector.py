"""按三个字符轻态的能量检验上下型零点选择，使用真实高度及原径向软化。"""
from functools import lru_cache
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx


def light_spectra(replay):
    shell = run_path(str(replay / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(replay / "ch39-native-light-projector-pump.py"))
    ctx.prec = 224
    blocks = list(pump["rotation_bases"](shell))
    adjoint = pump["adjoint"]

    @lru_cache(maxsize=None)
    def light(charge, index):
        matrix = shell["laplacian"](charge, acb.zeta_zero(index).imag)
        energies = []
        for block in blocks:
            spectrum = (adjoint(block)*matrix*block).eig(algorithm="rump")
            assert all(z.is_finite() and z.imag.contains(0) for z in spectrum)
            levels = sorted((z.real for z in spectrum), key=lambda e: float(e.mid()))
            assert levels[0] > 0 and all(levels[0] < e for e in levels[1:])
            energies.append(levels[0])
        return energies

    return light


def chebyshev_moments(energies, degree):
    previous = [arb(1) for _ in energies]
    x = [2*e/9-1 for e in energies]
    current = x
    values = []
    for _ in range(degree):
        values.append(sum(current))
        previous,current = current,[2*t*b-a for t,a,b in zip(x,previous,current,strict=True)]
    return values


def interior_radius(displacements):
    degree = len(displacements)-1
    matrix = arb_mat([[v[i] for v in displacements] for i in range(degree)]
                     + [[1]*(degree+1)])
    inverse = matrix.inv()
    weights = [inverse[i,degree] for i in range(degree+1)]
    assert all(w > 0 for w in weights)
    # 重心坐标在对面为零，用其原点值除以线性部分的范数取面距。
    return min((w/sum(inverse[i,j]**2 for j in range(degree)).sqrt()).lower()
               for i,w in enumerate(weights))


def verify():
    light = light_spectra(Path(__file__).resolve().parent)
    centers = {-4:2345666, 2:964890}
    for charge,first,runner in ((-4,73,-38),(2,-71,39)):
        center = centers[charge]
        def score(offset):
            return sum(e**2 for e in light(charge,center+offset))
        gaps = {i:score(0)-score(i) for i in range(-abs(first),abs(first)+1) if i}
        assert gaps[first] < 0 and all(g>0 for i,g in gaps.items() if i!=first)
        assert all(gaps[runner]<g for i,g in gaps.items() if 0<abs(i)<=64 and i!=runner)
        print("平方和",charge,"±64 次大值的偏移和差距",runner,gaps[runner].str(13),
              "首次超越的偏移和差距",first,(-gaps[first]).str(13))

    displacements = []
    for charge,index in ((-4,294),(-4,1411),(2,538),(2,1556)):
        target = chebyshev_moments(light(charge,centers[charge]),3)
        competitor = chebyshev_moments(light(charge,index),3)
        displacements.append([b-a for a,b in zip(target,competitor,strict=True)])
    radius = interior_radius(displacements)
    assert radius > arb(".00000753217")
    print("三次共同作用的矩位移凸包内半径",radius.str(14))


if __name__ == "__main__":
    verify()
