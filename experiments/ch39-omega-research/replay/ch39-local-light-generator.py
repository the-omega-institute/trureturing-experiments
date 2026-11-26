"""十九站点局域生成元实现无轻重跃迁时的最小重块及纯轻作用的距离下界。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


mix = run_path(str(Path(__file__).with_name("ch39-two-boundary-site-mixing.py")))
adjoint = mix["base"]["adjoint"]
sites = mix["shell"]["geo"]["sites"]


def identity(n):
    return acb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def forbidden(cutoff):
    return [(i,j) for i in range(19) for j in range(19)
            if sum((x-y)**2 for x,y in zip(sites[i], sites[j])) > cutoff]


def restriction(frame, edges):
    return acb_mat([[frame[i,a]*frame[j,b].conjugate() for a,b in product(range(3), repeat=2)]
                    for i,j in edges])


def positive(matrix):
    # 理论矩阵厄米，Sylvester 的顺序主子式判据。
    for n in range(1, matrix.nrows()+1):
        det = acb_mat([[matrix[i,j] for j in range(n)] for i in range(n)]).det()
        assert det.imag.contains(0) and det.real > 0


def heavy_gap(number):
    height = acb.zeta_zero(number).imag
    shell = mix["shell"]
    neutral = shell["laplacian"](0, height)
    center = sum(neutral[i,i].real for i in range(19)) / 18
    eigen = shell["laplacian"](1, height).eig(multiple=True)
    masses = [abs(z.real*(z.real-center)) for z in eigen]
    heavy_masses = sorted((x for x in masses if x > 9), key=float)
    assert len(heavy_masses) == 16
    gap = min((b-a).lower() for a,b in zip(heavy_masses, heavy_masses[1:]))
    assert gap > 0
    return gap


def verify(number):
    frame = mix["light_modes"](number)
    heavy = identity(19) - frame*adjoint(frame)
    edges = forbidden(804)
    assert len(edges) == 72
    actual = {(i,j) for i,j,_,_ in mix["shell"]["geo"]["edges"]}
    assert all(((i,j) in actual or (j,i) in actual) == ((i,j) not in edges)
               for i in range(19) for j in range(i))
    t = restriction(frame, edges)
    gram = acb_mat([[heavy[i,k]*heavy[l,j] for k,l in edges] for i,j in edges])
    # gram = R R* 半正定，区间解认证非奇异。
    inverse_t = gram.solve(t)
    optimal = adjoint(t)*inverse_t
    positive(optimal - identity(9)*arb(".00673")**2)
    positive(identity(9)*arb(".02161")**2 - optimal)
    gap = heavy_gap(number)
    print(number, "重质量的最小间隔及质量对易子下界系数", gap.str(9), (arb(".00673")*gap).str(9))

    for cutoff in (1152,1376):
        t = restriction(frame, forbidden(cutoff))
        # 前者有九个独立列，后者有六个独立行。
        matrix = adjoint(t)*t if cutoff == 1152 else t*adjoint(t)
        assert matrix.nrows() == (9 if cutoff == 1152 else 6)
        positive(matrix)

    if number == 28325:
        frames = mix["boundary_frames"](frame)
        controls, _ = mix["calibration"](frames)
        sum_bound = arb(0)
        for i,j,site,strength in controls:
            entry = strength*frames[site][i].conjugate()*frames[site][j]
            vec = acb_mat(9,1)
            vec[3*i+j,0] = entry
            vec[3*j+i,0] = -entry.conjugate()
            norm = (adjoint(vec)*optimal*vec)[0,0].real.sqrt()
            sum_bound += 2*norm
        print("五次局域 Cayley 作用造成的重空间酉矩阵偏移上界", sum_bound.str(12))
        assert sum_bound < arb(".006338")


if __name__ == "__main__":
    ctx.prec = 384
    for node in (28325,48067,58214):
        verify(node)
