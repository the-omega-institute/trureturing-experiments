"""原生十八个高模在径向边刚度 (d/12)^p、-20<=p<=20 时的一圈排除。

令 P=C+J/105。对共同质量尺，单个正规化的有限减法响应 F(m)
关于 log m 单调递减、斜率绝对值不超过一；这可逐项微分其核
log[(m²+z)/(m²+z0)]/2。将两个塔分别按质量递增配对，
F_b-F_a <= B_ab=sum_j max(0, log(m_aj/m_bj))。
同谱时标量的角核分布 y³ <= (8y-y³)/7，故 F_H<=F_V。
三个伴随与 Higgs 在 P 中的贡献至多 B_(0,3)/(5 pi)。
三代 Weyl 系数按正、负部分配对，用相同的 B_ab 界。

对增量 h>=0，各边都是半正定，L_q(p)<=L_q(p+h)<=67^h L_q(p)。
因此每个有序 log m 在 [0,h log(67)/2] 内变化，每对谱的 B_ab
在每一侧至多变化 9 h log(67)。下面的 161 点覆盖整个闭区间。
所用轻场运行来自中心输入，物理匹配窗口为 10^13 至 10^17 GeV。
"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx, fmpz_mat

census = run_path(str(Path(__file__).with_name("ch39-shell-charged-channel-census.py")))
soft = census["soft"]
charges = (0, 1, 4, 2, 3, 6)


def distance_blocks():
    h = acb.zeta_zero(11431).imag
    geo = soft["geo"]
    phases = geo["phases"](h)
    soft["star"]["charged"]["geometry"]()  # integer proof of cyclic invariance
    bases = [census["rotation_basis"](s) for s in range(3)]
    matrices = {}
    for i, j, d, z in geo["edges"]:
        if d not in matrices:
            matrices[d] = [acb_mat(19, 19) for _ in charges]
        for q, matrix in zip(charges, matrices[d]):
            link = acb(1) if soft["pole"] in (i, j) else (z*phases.get(d, 1))**q
            matrix[i, i] += 1
            matrix[j, j] += 1
            matrix[i, j] = -link
            matrix[j, i] = -link.conjugate()
    assert min(matrices) == 12 and max(matrices) == 804
    blocks = {d:[[basis.conjugate().transpose()*matrix*basis for basis in bases]
                 for matrix in row] for d, row in matrices.items()}
    base, neutral, _ = soft["vector_pencil"](h)
    flat = base+neutral
    integer = fmpz_mat([[int(flat[i,j].real.unique_fmpz()) for j in range(19)] for i in range(19)])
    unweighted = [z.real for z, mult in integer.charpoly().complex_roots() for _ in range(mult)]
    return blocks, unweighted


def log_masses(blocks, unweighted, p):
    scales = {d:(arb(d)/12)**p for d in blocks}
    result = []
    for idx, q in enumerate(charges):
        if p == 0 and q == 0:
            values = unweighted
        else:
            values = []
            for sector, size in enumerate((7, 6, 6)):
                block = sum((scales[d]*row[idx][sector] for d, row in blocks.items()),
                            acb_mat(size, size))
                values.extend(z.real for z in block.eig(multiple=True))
        # Separate sorting of endpoints encloses order statistics even for repeated roots.
        lower = sorted(v.lower() for v in values)
        upper = sorted(v.upper() for v in values)
        assert len(lower) == 19 and upper[0] < lower[1] and lower[1] > 0
        result.append([lower[k].union(upper[k]).log()/2 for k in range(1, 19)])
    return result


def spread(first, second):
    return sum((u-v).max(arb(0)) for u, v in zip(first, second))


def verify():
    ctx.prec = 330
    blocks, unweighted = distance_blocks()
    pi = arb.pi()
    # Single-generation beta in P, ordered by |6Y|=(1,4,2,3,6).
    weights = [arb(i)/105 for i in (-172,128,84,-84,44)]
    assert sum(weights).contains(0)
    total = sum(w for w in weights if w > 0)
    em, weak = arb("127.930"), arb(".23129")
    y, w, c = 3*(1-weak)*em/5, weak*em, 1/arb(".1180")
    log_min = (arb("1e13")/arb("91.1876")).log()
    baseline = y-3*w+2*c+(5*y-3*w-2*c)/105-log_min/(105*pi)
    assert baseline < arb("-11.08713")
    peak = arb(0)
    for index in range(-80, 81):
        spectra = log_masses(blocks, unweighted, arb(index)/4)
        scalar = spread(spectra[0], spectra[4])/(5*pi)
        weyl = sum(wi*(-wj)*spread(spectra[i+1], spectra[j+1])
                   for i, wi in enumerate(weights) for j, wj in enumerate(weights)
                   if wi > 0 and wj < 0)
        peak = peak.max(scalar+3*weyl/(2*pi*total))
    cell_error = (1/(5*pi)+3*total/(2*pi))*9*arb(67).log()/8
    bound = baseline+peak+cell_error
    assert peak < arb("1.095135") and cell_error < arb("5.807680")
    assert bound < arb("-4.18431") and -bound*35/106 > arb("1.381")
    print("161 radial exponents, p=-20 through 20; maximum grid correction", peak.upper().str(13))
    print("Largest SM projection and cell allowance", baseline.upper().str(13), cell_error.upper().str(13))
    print("Uniform P upper and minimum inverse-coupling spread", bound.upper().str(13),
          (-bound*35/106).lower().str(13))


if __name__ == "__main__":
    verify()
