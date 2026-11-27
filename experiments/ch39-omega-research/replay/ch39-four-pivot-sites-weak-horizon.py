"""只在极点与三条枢轴射线站位施酉，从第七壳起精确消去弱差的最大长度。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(matrix):
    return matrix.conjugate().transpose()


def norm(matrix):
    return sum(abs(z)**2 for z in matrix.entries()).sqrt()


def column(matrix, index):
    return acb_mat([[row[index]] for row in matrix.tolist()])


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    sieve = lib("six-cyclic-line-sieve")
    candidates = lib("seven-cyclic-weak-horizon")["CANDIDATES"][:3]
    # 沿用逐壳弱核界：若每条候选线都至少 3000 GeV，则该壳弱差严格为负。
    _, (_, _, model) = sieve["prepare_shells"](root, fifth, columns, levels, candidates, arb(0))
    forms, _ = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    sites = [graph["soft"]["pole"], *graph["boundary_sites"]()]
    assert sites == [14, 8, 7, 6]
    gram = native["gram_cholesky"](root)
    floor, heavy, top = levels
    ctx.prec = 256
    planck, mz = arb("1.22089e19"), arb("91.1876")
    source = fifth * columns * acb_mat([[(heavy-floor).sqrt(), 0], [0, (top-floor).sqrt()]])
    local = acb_mat([source.tolist()[i] for i in sites])
    spectra = []
    for j, ids in enumerate(candidates):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, 2_306_888_100_092+j)
        rows = adj(basis)
        restricted = acb_mat([[row[i] for i in sites] for row in rows.tolist()])
        old = rows * source
        weights = [z.real for z in native["readout"](basis, range(19), sites, gram[3], 3)]
        assert all(w.overlaps(v) for w, v in zip(weights, model[j]["w"], strict=True))
        # 原生谱简单且图与三循环对易，平凡字符分量非零即标识整个谱行属于该字符。
        assert all(norm(acb_mat([rows.tolist()[i]]) * forms[0]) > arb(".9") for i in ids)
        spectra.append((restricted, old, weights))

    # 从每壳选一条候选线，在其四站限制上，三个射线分量严格相同。
    bounds = []
    for ids in product(*candidates):
        pole, ray, before = [], [], []
        for (restricted, old, _), i in zip(spectra, ids, strict=True):
            pole.append(restricted[i, 0])
            ray.append(sum(restricted[i, k] for k in (1, 2, 3))/3)
            before.append(acb_mat([old.tolist()[i]]))
        minors = [pole[(j+1)%3]*ray[(j+2)%3] - pole[(j+2)%3]*ray[(j+1)%3]
                  for j in range(3)]
        total = sum(map(abs, minors))
        assert total > 0
        preserved = minors[0]*before[0] + minors[1]*before[1] + minors[2]*before[2]
        bounds.append(((floor + (norm(preserved)/total)**2)*planck).lower())
    assert len(bounds) == 18 and min(bounds) > arb("3.25e14")
    minimum_row = (min(bounds)/planck - floor).sqrt()
    low_row = (arb(3000)/planck-floor).sqrt()
    separation = (minimum_row-low_row)/(top-floor).sqrt()
    assert separation > arb(".0059357")

    # 第七、八壳分别用第十、十六线：固定其激发振幅，并在四站内补齐源 Gram。
    selected = (9, 15)
    constraints = acb_mat([s[0].tolist()[i] for s, i in zip(spectra[:2], selected, strict=True)])
    before = acb_mat([s[1].tolist()[i] for s, i in zip(spectra[:2], selected, strict=True)])
    outside = before - constraints*local
    dual = adj(constraints)*(constraints*adj(constraints)).inv()
    eye = acb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    kernel = eye - dual*constraints
    u = column(kernel, 2)
    u = u/norm(u)
    v = column(kernel, 3)
    v = v - u*(adj(u)*v)[0, 0]
    v = v/norm(v)
    null_frame = acb_mat([x+y for x, y in zip(u.tolist(), v.tolist(), strict=True)])
    endpoints = [(arb("178.45518"), arb("178.45519")), (arb("895.56334"), arb("895.56335"))]
    target = acb_mat([[((lo.union(hi))/planck-floor).sqrt(), 0] for lo, hi in endpoints])
    lift = dual*(target-outside)
    remaining = adj(local)*local-adj(lift)*lift
    assert remaining[0, 0].real > 0
    a = remaining[0, 0].real.sqrt()
    b = remaining[0, 1]/a
    c2 = remaining[1, 1].real-abs(b)**2
    assert c2 > 0
    image = lift + null_frame*acb_mat([[a, b], [0, c2.sqrt()]])
    weak = native["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    reference = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    p0, p1 = (mz/planck)**2, (arb("1.176397e14")/planck)**2
    for (restricted, old, weights), i, ends in zip(spectra[:2], selected, endpoints, strict=True):
        changed = old + restricted*(image-local)
        minimum = min((floor+sum(abs(z)**2 for z in row)).lower()
                      for j, row in enumerate(changed.tolist()) if j != i)
        assert minimum > arb("1e-5")
        tail = (large(minimum**2, p0, p1)-reference).abs_upper()
        tail *= sum(abs(w).upper() for j, w in enumerate(weights) if j != i)
        faces = [reference + weights[i]*(weak((q/mz)**2)-reference) + arb(0, tail)
                 for q in ends]
        assert faces[0] > 0 > faces[1]
        print("低线与两弱差面的包围", i+1, *[x.str(12) for x in faces], flush=True)
    print("三壳候选组合数与至少一条所选线的截止下界（GeV）",
          len(bounds), min(bounds).str(14), flush=True)


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root/"ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))
