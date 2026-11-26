"""用一个外站的谱极点检验三枚算术枢轴的共同质量标尺。"""
from fractions import Fraction
from functools import cache
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


@cache
def phase(n, index):
    return acb(0, arb.pi()*(2*index+1)/n).exp()


@cache
def variation(n):
    return 4*arb(2).sqrt()*(arb.pi()/(2*n)).sin() + 6*arb("1e-8")


def adj(matrix):
    return matrix.conjugate().transpose()


def verify(root):
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    native = run_path(str(root/"ch39-native-light-projector-pump.py"))
    levels = run_path(str(root/"ch39-three-site-readout.py"))
    pivots = run_path(str(root/"ch39-finite-rh-pivot-mass.py"))
    ctx.prec = 128
    small_minor = pivots["determinant"](range(2))
    medium_minor = pivots["determinant"](range(3))
    large_minor = pivots["determinant"](range(4))
    ratio_squared_lower = medium_minor[0]**2/(small_minor[1]*large_minor[1])
    assert min(small_minor[0],medium_minor[0],large_minor[0]) > 0
    assert ratio_squared_lower > Fraction(49,5)**2
    bases = list(native["rotation_bases"](shell))
    for q in (1,-4,2,-3,6):
        base, clocks = native["clock_matrices"](shell,q)
        blocks = [[adj(b)*m*b for m in (base,*clocks)] for b in bases]
        for b in bases:
            assert [{i for i in range(19) if not b[i,j].is_zero()} for j in (0,1)] == [
                {0,3,4}, {1,2,5}]
        leaves = finest = 0
        residue_bound, ratio_bound = arb(1),arb(0)
        todo = [(8,i,j) for i in range(8) for j in range(8)]
        while todo:
            n,i,j = todo.pop()
            perturbation = variation(n)
            energies, next_levels, weights = [], [], []
            for constant,first,second in blocks:
                off = phase(n,i)*first + phase(n,j)*second
                v,energy,gap = levels["lowest"](constant+off+adj(off))
                if not gap>perturbation:
                    break
                energies.append(energy)
                next_levels.append((energy+gap-perturbation).lower())
                amplitude = min(abs(v[t]).lower() for t in (0,1))-perturbation/(gap-perturbation)
                if not amplitude>0:
                    break
                weights.append((amplitude**2/3).lower())
            if len(weights)==3:
                denominator = sorted(z.lower() for z in energies)[1]-perturbation
                upper = max(z.upper() for z in energies)+perturbation
                if (denominator>0 and min(weights)>arb(1)/48 and min(next_levels)>upper
                        and upper<arb(6)/5*denominator):
                    leaves += 1
                    finest = max(finest,n)
                    residue_bound = min(residue_bound, min(weights))
                    ratio_bound = max(ratio_bound,(upper/denominator).upper())
                    continue
            assert n<4096, (q,n,i,j)
            todo.extend((2*n,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
        print("荷、角方格数、最细每圈格数、单站最低谱权下界、第三与第二极点比上界",
              q,leaves,finest,residue_bound.str(10),ratio_bound.str(10),flush=True)
    assert Fraction(49,5)/Fraction(6,5)>8
    print("有理 Gram 行列式包围给出相邻枢轴比大于 9.8")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
