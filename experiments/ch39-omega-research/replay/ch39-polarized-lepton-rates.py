"""由极化带电轻子的类空约化散射率反演中性交换，核查符号歧义。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, fmpq_mat

light = run_path(str(Path(__file__).with_name("ch39-light-weyl-current-basis.py")))
poles, q = light["poles"], light["q"]
entries = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
channels = ("AA", "AB", "BB", "AC", "AD", "BC", "BD", "CD")
basis_names = ("AA", "AB", "BB", "AC", "BD", "CD")
rate_errors = ("1e-14", "1e-9", "1e-9")


def charge(name, x, conjugation):
    d, r, sign = light["charges"][name]
    return tuple(conjugation*v for v in (d, r, sign*x))


leptons = dict(A=charge("L", 1, 1), B=charge("ec", 1, -1),
               C=charge("L", -2, 1), D=charge("ec", -2, -1))


def rows(names):
    result = []
    for left, right in names:
        v, w = leptons[left], leptons[right]
        result.append([v[i]*w[j] if i == j else v[i]*w[j]+v[j]*w[i] for i, j in entries])
    return fmpq_mat(result)


def as_matrix(vector):
    a, b, c, d, e, f = vector.entries()
    return arb_mat([[a, d, e], [d, b, f], [e, f, c]])


def column(g):
    return arb_mat([[g[i, j]] for i, j in entries])


def signs(names):
    return product(*[(1,) if a == b else (1, -1) for a, b in names])


def has_negative_minor(g):
    if any(g[i, i] < 0 for i in range(3)):
        return True
    if any(g[i, i]*g[j, j] - g[i, j]**2 < 0 for i, j in combinations(range(3), 2)):
        return True
    return g.det() < 0


def sign_certificate(graph, error):
    basis = rows(basis_names)
    inverse = arb_mat(basis.inv())
    central = (arb_mat(basis)*column(graph)).entries()
    extra = arb_mat(rows(("AD",)))
    extra_rate = (extra*column(graph))[0, 0]**2
    assert extra_rate > 0
    relative_window = 1+arb(0, error.upper())
    extra_window = extra_rate*relative_window
    moduli = [(a**2*relative_window).sqrt() for a in central]
    survivors = []
    for sign in signs(basis_names):
        signed = arb_mat([[s*m] for s, m in zip(sign, moduli, strict=True)])
        candidate = inverse*signed
        squared = (extra*candidate)[0, 0]**2
        if squared.overlaps(extra_window) and not has_negative_minor(as_matrix(candidate)):
            survivors.append(sign)
    assert all(abs(a) > 0 for a in central)
    reference = tuple(1 if a > 0 else -1 for a in central)
    assert survivors == [reference]
    return reference


def six_rate_ambiguities(graphs):
    deficient, ambiguous = 0, 0
    for names in combinations(channels, 6):
        basis = rows(names)
        if basis.rank() < 6:
            deficient += 1
            continue
        inverse = arb_mat(basis.inv())
        found = False
        for graph in graphs:
            exact = (arb_mat(basis)*column(graph)).entries()
            for sign in signs(names):
                if all(s*a > 0 for s, a in zip(sign, exact, strict=True)):
                    continue
                alternate = as_matrix(inverse*arb_mat([[s*abs(a)] for s, a in zip(sign, exact, strict=True)]))
                if light["mat"]["positive"](alternate) and any(s*a < 0 for s, a in zip(sign, exact, strict=True)):
                    found = True
                    break
            if found:
                break
        assert found, names
        ambiguous += 1
    assert (deficient, ambiguous) == (13, 15)
    return deficient, ambiguous


def verify():
    for i in range(3):
        assert leptons["A"][i]+leptons["D"][i] == leptons["B"][i]+leptons["C"][i]
    unseen = fmpq_mat([[1], [0], [9], [0], [3], [0]])
    assert rows(channels[:-1])*unseen == fmpq_mat(7, 1)
    assert (rows(("CD",))*unseen)[0, 0] == 81
    assert rows(channels[:-1]).rank() == 5

    time = poles["study"]["nominal"][2]
    models = [light["graphs"](logs, time) for logs in poles["retune"]["log_sets"]]
    assert len(models) == 3
    assert all(light["mat"]["positive"](g) for model in models for g in model)
    print("six-rate subsets, deficient and ambiguous", six_rate_ambiguities(
        [g for model in models for g in model]), flush=True)

    for case, model in enumerate(models):
        recovered = [sign_certificate(g, q(err)) for g, err in zip(model, rate_errors, strict=True)]
        assert [s[basis_names.index("AC")] for s in recovered] == [1, -1, -1]
        basis = rows(basis_names)
        limits = [light["thresholds"](model, basis, n) for n in (9, 10)]
        rooms = [min(column) for column in zip(*limits, strict=True)]
        print("case", case, "signs by momentum", recovered, "amplitude tolerances",
              *(v.str(10) for v in rooms), flush=True)
        for room, err in zip(rooms, rate_errors, strict=True):
            delta = q(err)
            assert delta/(1+(1-delta).sqrt()) < room


if __name__ == "__main__":
    verify()
