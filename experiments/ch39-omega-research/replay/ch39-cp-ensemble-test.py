"""已校准的 CKM 概率之外，检验按七个实际背景平均的 CP 奇不变量。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, arb_poly, arb_series, ctx


root = Path(__file__).parent
two = run_path(str(root / "ch39-two-ckm-probabilities.py"))
flavor = two["invariant"]["flavor"]


def cp_data(data, matrix):
    e1, e2, e3 = data[:3]
    discriminant = e1**2*e2**2-4*e2**3-4*e1**3*e3-27*e3**2+18*e1*e2*e3
    unit = arb(".01401")**2
    nodes = [arb(j)/4 for j in (-3, -1, 1, 3)]
    values = []
    for c in nodes:
        y = matrix(two["invariant"]["phase"](c))
        h = (y*y.conjugate().transpose())*(1/unit)
        numerator = (h[0,1]*h[1,2]*h[2,0]).imag
        values.append([numerator/(1-c*c).sqrt()])
    fit = arb_mat([[c**i for i in range(4)] for c in nodes]).solve(arb_mat(values))
    return arb_poly(fit.entries()), discriminant


def cp_value(c, polys, jet=False):
    x = arb_series([c, 1], prec=3) if jet else c
    radicand = 1-x*x
    sine = radicand.sqrt() if jet else radicand.nonnegative_part().sqrt()
    numerator, discriminant = [two["evaluate"](p, x) for p in polys]
    return sine*numerator/discriminant.sqrt()


def certify_supporting_line(data, polys, reference):
    up, _ = two["probabilities"](reference, data)
    cp = cp_value(reference, polys, jet=True)
    slope = cp[1]/up[1]
    margin = arb(".0006")
    interior = reference+arb(0, arb(".0029"))
    for j in range(32):
        ends = [reference + arb(".003")*(arb(i)/16-1) for i in (j, j+1)]
        subdomain = ends[0].union(ends[1])
        nearby_up, _ = two["probabilities"](subdomain, data)
        assert slope*nearby_up[2]-cp_value(subdomain, polys, jet=True)[2] > margin
    pending = [(arb(-1), arb(1))]
    while pending:
        a, b = pending.pop()
        domain = a.union(b)
        if interior.contains(domain):
            continue
        local_up, _ = two["probabilities"](domain, data)
        gap = cp[0]+slope*(local_up[0]-up[0])-cp_value(domain, polys)
        if gap - margin*(domain-reference)**2 > 0:
            continue
        mid = (a+b)/2
        assert b-a > arb("1e-10")
        pending.extend(((a,mid), (mid,b)))
    return slope


def matrix_readings(matrix, phase):
    y = matrix(phase)
    h = y*y.conjugate().transpose()
    eigen = sorted(v.real for v in h.eig(algorithm="rump"))
    assert 0 < eigen[0] < eigen[1] < eigen[2]
    eye = acb_mat([[int(i==j) for j in range(3)] for i in range(3)])
    readings = []
    for j, mass in enumerate(eigen):
        projector = eye
        for k, other in enumerate(eigen):
            if j != k:
                projector = projector*(h-other*eye)/(mass-other)
        readings.extend(projector[i,i].real for i in range(3))
    cp = (h[0,1]*h[1,2]*h[2,0]).imag/((eigen[1]-eigen[0])*(eigen[2]-eigen[0])*(eigen[2]-eigen[1]))
    return readings, eigen, cp


def verify():
    ctx.prec, ctx.cap = 384, 3
    data = two["spectral_data"]()
    matrix = flavor["reference_matrix"]()
    phase, c0 = flavor["singlet"](11661)
    ref, masses, cp0 = matrix_readings(matrix, phase)
    polys = cp_data(data, matrix)
    assert cp_value(c0, polys).overlaps(cp0)
    slope = certify_supporting_line(data, polys, c0)
    joint = two["joint"]
    model = joint["model"]
    times = model["shared"]["clock"]["clock"]["times"]()
    _, _, lambdas, slopes = model["shared"]["flavor_data"]()
    columns, readings = [], []
    for n in joint["indices"]:
        phase, coordinate = flavor["singlet"](n)
        p, m, cp = matrix_readings(matrix, phase)
        assert cp_value(coordinate, polys).overlaps(cp)
        readings.append((p, m, cp))
        x = (joint["scale"]/acb.zeta_zero(n).imag)**2
        columns.append([1, *[1/(1+j*j*x) for j in range(1,5)], model["gradient"](x,lambdas,slopes), p[6]])
    weights = arb_mat(list(zip(*columns, strict=True))).solve(arb_mat([[v] for v in [*times[:5],0,ref[6]]])).entries()
    assert all(w > 0 for w in weights)
    differences = [sum(w*v[0][j] for w, v in zip(weights, readings))-target for j, target in enumerate(ref)]
    worst = max(abs(v).upper() for v in differences)
    mass_error = max(abs(sum(w*v[1][j] for w,v in zip(weights,readings))/target-1).upper()
                     for j,target in enumerate(masses))
    cp_deficit = 1 - sum(w*v[2] for w, v in zip(weights, readings))/cp0
    assert worst < arb("9.214e-9") and mass_error < arb("2.683e-10")
    assert arb(".01409299") < cp_deficit < arb(".01409300")
    print("CP 切线系数与七点平均 CP 相对亏损", slope.str(11), cp_deficit.str(12))
    print("九个跃迁概率的最大绝对误差上界、三个质量平方的最大相对误差上界", worst.str(12), mass_error.str(12))


if __name__ == "__main__":
    verify()
