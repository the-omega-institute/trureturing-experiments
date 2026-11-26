"""用已有轻 Weyl 场的自交换重建有理动量协议所需的读数。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, fmpq, fmpq_mat

scheme = run_path(str(Path(__file__).with_name("ch39-rational-momentum-calibration.py")))
cal, mat, poles, q = (scheme[key] for key in ("cal", "mat", "poles", "q"))
charges = poles["study"]["retuned"]["current"]["weyl_charges"]
family, higgs, hyper = (0, 0, 1), (3, -1, 2), (1, 3, 0)


def lift(vector):
    a, b, c = vector
    return [a*a, b*b, c*c, 2*a*b, 2*a*c, 2*b*c]


def current_bases():
    tagged, inclusive = {}, {}
    for name, (d, r, sign) in charges.items():
        for x in (1, -2):
            tagged[name, x] = lift((d, r, sign*x))
        inclusive[name] = [fmpq(2*a+b, 3) for a, b in
                           zip(tagged[name, 1], tagged[name, -2], strict=True)]

    leptons = [row for (name, _), row in tagged.items() if name in ("L", "ec")]
    lepton_matrix = fmpq_mat(leptons)
    assert lepton_matrix.rank() == 4
    unseen = fmpq_mat([[4], [0], [-18], [0], [3], [0]])
    assert lepton_matrix*unseen == fmpq_mat(4, 1)
    assert (fmpq_mat([lift(family)])*unseen)[0, 0] == -18

    quarks = [key for key in tagged if key[0] not in ("L", "ec")]
    options = {}
    for pair in combinations(quarks, 2):
        matrix = fmpq_mat(leptons + [tagged[key] for key in pair])
        if matrix.rank() == 6:
            options[pair] = matrix
    assert len(options) == 14
    assert (("Q", 1), ("Q", -2)) not in options

    assert fmpq_mat(list(inclusive.values())).rank() == 4
    assert all(row[-2:] == [0, 0] for row in inclusive.values())
    averaged = [inclusive[name] for name in ("Q", "L", "uc", "dc")]
    two_tags = fmpq_mat(averaged + [tagged["uc", -2], tagged["dc", 1]])
    assert two_tags.rank() == 6
    return options[(("uc", -2), ("dc", 1))], two_tags


def graphs(logs, time):
    (length, z, rho, eta, c, masses), _, _ = poles["retune"]["predictions"](logs, time)
    e = scheme["anchor"]/length
    h = mat["flow"]["running"]["hypercharge_window"](time + scheme["anchor"].log()/2).mid()**2

    def bare(exchange):
        return exchange*(1/rho) + mat["hypercharge"]*h

    result = [bare(mat["massive_exchange"](e, z, c/rho, eta/rho))]
    for log_x in scheme["logs_x"]:
        exchange, polarization, _ = cal["shape"]["probe"](log_x, e, z, c/rho, eta/rho, masses[0])
        tree = bare(exchange)
        result.append((cal["eye"] - tree*polarization).inv()*tree)
    return result


def thresholds(graph, basis, n):
    inverse = basis.inv()
    entries = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)]
    read = [(arb_mat(basis)*arb_mat([[g[i, j]] for i, j in entries])).entries() for g in graph]
    assert all(abs(y) > 0 for point in read for y in point)

    def nominal(point, vector):
        col = arb_mat([[v] for v in vector])
        return cal["form"](col, graph[point], col)

    def coeff(vector):
        return [arb(a) for a in (fmpq_mat([lift(vector)])*inverse).entries()]

    def cost(point, row):
        return sum(abs(a*y) for a, y in zip(row, read[point], strict=True))

    def same(vector, eps):
        numerator, denominator = nominal(0, vector), nominal(0, family)
        ratio = numerator/denominator
        a, b = coeff(vector), coeff(family)
        caps = [denominator/cost(0, b)]
        for sign in (-1, 1):
            residual = [p-ratio*(1+sign*eps)*h for p, h in zip(a, b, strict=True)]
            caps.append(numerator*eps/cost(0, residual))
        return min(cap.lower() for cap in caps)

    # 三个动量共用归一标尺。分子和分母各分配一半跨动量误差预算。
    def across(vector, point, eps):
        row = coeff(vector)
        cap0 = eps*nominal(0, vector)/(2*(1+eps)*cost(0, row))
        cap1 = eps*nominal(point, vector)/(2*cost(point, row))
        return cap0.lower(), cap1.lower()

    shape_room = same(higgs, q("1e-10"))
    hyper_room = same(hyper, q("1e-6"))
    f1, f2 = [across(family, point, q("2e-8")) for point in (1, 2)]
    mix = across((3*n+1, -n, 2*n), 2, q("2e-5"))
    return [min(shape_room, hyper_room, f1[0], f2[0], mix[0]), f1[1], min(f2[1], mix[1])]


def verify():
    tagged, two_tags = current_bases()
    time = poles["study"]["nominal"][2]
    models = [graphs(logs, time) for logs in poles["retune"]["log_sets"]]
    assert len(models) == 3
    cases = (("six tagged", tagged, ("5e-14", "2e-9", "2e-9")),
             ("four averages and two tags", two_tags, ("2e-14", "5e-9", "5e-9")))
    for name, basis, precision in cases:
        limits = [thresholds(graph, basis, n) for graph in models for n in (9, 10)]
        bounds = [min(column) for column in zip(*limits, strict=True)]
        assert all(b > q(e) for b, e in zip(bounds, precision, strict=True))
        print(name, "det", basis.det(), "relative errors by momentum",
              *(v.str(12) for v in bounds), flush=True)


if __name__ == "__main__":
    verify()
