"""在三个整圆与三个二维环面上排除至多两条最短站位边的实旋转。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb


def experiment(root, fifth, columns, levels):
    shared = run_path(str(root / "ch39-fifth-seventh-shortest-site-gates.py"))
    baseline, cosine, sine, _, _, weights, heat, field = shared["geometry"](root, fifth, columns, levels)
    low, heavy, highest = levels
    increases = (heavy-low, highest-low)
    coefficient = -sum(g for g, _, _ in field)/(8*arb.pi())
    for count in (1, 2):
        for axes in combinations(range(3), count):
            fixed = baseline
            for edge in axes:
                fixed = fixed-cosine[edge]

            # 全角域的逐谱线振幅上界。热核质量损失随平方截止单调增加。
            errors = []
            for i, weight in enumerate(weights):
                caps = [abs(fixed[i,j])+sum(abs(cosine[e][i,j])+abs(sine[e][i,j]) for e in axes)
                        for j in range(2)]
                amplitude = min(highest.upper(), (low+sum(d*c*c for d,c in zip(increases,caps,strict=True))).upper())
                x = amplitude**2
                loss = x*x*sum(abs(g)*(arb(1)/2-heat(3,z/x)) for g,_,z in field)/(4*arb.pi())
                errors.append(abs(weight)*loss)
            error = sum(errors).upper()
            bound = arb(-1) if count == 1 else arb("-.49")
            cells = [tuple((arb(0),2*arb.pi()) for _ in axes)]
            leaves = 0
            ceiling = None
            while cells:
                cell = cells.pop()
                overlaps = fixed
                for e, (left,right) in zip(axes,cell,strict=True):
                    s,c = left.union(right).sin_cos()
                    overlaps = overlaps+c*cosine[e]+s*sine[e]
                amplitudes = [low+sum(d*abs(overlaps[i,j])**2 for j,d in enumerate(increases))
                              for i in range(19)]
                vacuum = coefficient*sum(c*a**4 for c,a in zip(weights,amplitudes,strict=True))
                upper = vacuum.upper()+error
                if upper < bound:
                    leaves += 1
                    ceiling = upper.upper() if ceiling is None else max(ceiling,upper.upper())
                    continue
                j = max(range(count),key=lambda k:(cell[k][1]-cell[k][0]).upper())
                left,right = cell[j]
                assert right-left > arb("1e-6")
                middle = (left+right)/2
                for part in ((left,middle),(middle,right)):
                    cells.append(cell[:j]+(part,)+cell[j+1:])
            print("活动边、闭盒数、质量修正绝对上界和真空项统一上界",
                  tuple(e+1 for e in axes),leaves,error.str(10),ceiling.str(13),flush=True)


def verify(root):
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
