"""固定源端雅可比与二阶核，用预登记的二十四枚目标重谱检验能标预报。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

ROOT = Path(__file__).resolve().parent
TARGETS = tuple(range(2346028, 2346075, 2))


def matrix_row(mat, index):
    return arb_mat([[mat[index, j] for j in range(mat.ncols())]])


def quadratic_part(moments):
    return arb_mat([[x] for x in [*moments.entries()[:5], moments.entries()[-1]]])


def high_window(physical, energy):
    if not physical:
        return False
    low, high = arb(10)**13, arb(10)**17
    assert not (energy-low).contains(0) and not (energy-high).contains(0)
    return low < energy < high


def predict():
    repair = run_path(str(ROOT / "ch39-one-share-gravity-repair.py"))
    load = repair["load"]
    _, kernel, moments, changes = repair["source_data"]()
    numbers = (repair["SOURCE"], *TARGETS)
    _, lib, spectra = load("two-zero-cutoff-gravity")["gravity_fields"](ROOT, numbers, tight_neutral=True)
    model, _, _, initial, _ = lib["chiral"]["source"]["mean"]["inputs"]()
    fields = [(0,0,0), *load("signed-gravity-unification")["gauge_weights"](model)]
    beta = [b/(2*arb.pi()) for b in model["beta"]]
    t = (beta[0]-beta[2])/(beta[0]-beta[1])
    directions = arb_mat([[1,0,0,0,0], [0,1,0,0,0], [0,0,1-t,t,-1]])
    src_rows = directions * kernel["series"](spectra[0], fields)
    deriv = []
    for k in range(3):
        row = matrix_row(src_rows, k)
        slope = (repair["dilate"](row, arb(1), True)*moments)[0,0]
        deriv.append([slope, *[(row*changes[line][0])[0,0] for line in (16,301)]])
    source_inverse = arb_mat(deriv).inv()
    source_coeff = kernel["series"](spectra[0], fields, order=2)
    m = quadratic_part(moments)
    m16, m301 = (quadratic_part(changes[i][0]) for i in (16,301))
    results = {}
    for n, rows in zip(TARGETS, spectra[1:], strict=True):
        coeff = kernel["series"](rows, fields, order=2)
        ds, d16, d301 = (-(source_inverse * directions * (coeff-source_coeff) * m)).entries()
        s = 1+ds
        assert s > 0
        target_moments = m+d16*m16+d301*m301
        start = []
        for i in range(3):
            transported = repair["dilate"](matrix_row(coeff, 2+i), s)
            start.append(initial[i]+(transported*target_moments)[0,0])
        ell = (start[0]-start[2])/(beta[0]-beta[2])
        energy = model["mz"] * ell.exp()
        endpoint = start[2]-beta[2]*ell
        p, q = changes[16][1]+d16, changes[301][1]+d301
        # 沿用已登记的一阶预报判据。在近似交点读取首尾两端与投影公共值。
        checks = [s, 1-arb(".007")*s, p, 1-p, q, 1-q, ell,
                  (arb("1.22089e19")/model["mz"]).log()-ell, *start, endpoint]
        assert all(not x.contains(0) for x in checks)
        physical = all(x > 0 for x in checks)
        results[n] = physical, energy
        print("二阶预报", n, physical, high_window(physical,energy), energy.str(12), flush=True)
    return results


def verify():
    forecasts = predict()
    certify = run_path(str(ROOT / "ch39-two-source-share-crossings.py"))["verify"]
    certified = certify(forward=True, targets_override=TARGETS)
    assert list(forecasts) == list(certified) == list(TARGETS)
    misclassified = {"物理": [[],[]], "高能": [[],[]]}
    physical_errors = []
    for n in TARGETS:
        proposed, estimate = forecasts[n]
        accepted, energy = certified[n]
        in_window = high_window(accepted,energy)
        for label, guess, actual in (("物理", proposed, accepted),
                                     ("高能", high_window(proposed,estimate), in_window)):
            if guess != actual:
                misclassified[label][int(actual)].append(n)
        print("对照", n, proposed, accepted, high_window(proposed,estimate), in_window,
              estimate.str(12), energy.str(12), flush=True)
        if accepted:
            physical_errors.append(abs(estimate/energy-1).upper())
    for label,(false_positives,false_negatives) in misclassified.items():
        print(label, "一致", len(TARGETS)-len(false_positives)-len(false_negatives),
              "假阳性", false_positives, "假阴性", false_negatives)
    assert [n for n,(accepted,_) in certified.items() if accepted] == [2346042,2346054]
    assert all(not e for pair in misclassified.values() for e in pair)
    assert not any(high_window(ok, energy) for ok, energy in certified.values())
    assert max(physical_errors) < arb(".052368")
    print("正例最大相对能标误差", max(physical_errors).str(12))


if __name__ == "__main__":
    verify()
