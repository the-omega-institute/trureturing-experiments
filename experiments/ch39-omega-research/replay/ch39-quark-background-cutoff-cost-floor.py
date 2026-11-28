"""在所有普朗克有界截止测度上给两个实际夸克背景的共同引力矩以负量下界。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify():
    root = Path(__file__).resolve().parent
    spectra = run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))
    _, lib, (spectrum, reference) = spectra["gravity_fields"](
        root, (27_804_036, 2_345_666), tight_neutral=True)
    heat = lib["expint"]
    pi = arb.pi()
    total_alpha = sum(a for _, a, _ in spectrum)
    delta_moment = sum(g*(s-t) for (g, _, s), (_, _, t) in zip(spectrum, reference, strict=True))
    assert not delta_moment.contains(0)
    factor = 2*total_alpha/delta_moment

    # H=K_s−factor·(W_s−W_t) 的线性主项抵消。每项均在 0 与给出的有符号端点间。
    endpoints = [-a*(1-heat(2, s))/(2*pi) for _, a, s in spectrum]
    for (g, a, s), (h, b, t) in zip(spectrum, reference, strict=True):
        assert g == h and (a-b).contains(0)
        endpoints.append(factor*g*(s-t)*(1-heat(2, s.max(t)))/(4*pi))
    lower = sum(v.min(0) for v in endpoints).lower()
    upper = sum(v.max(0) for v in endpoints).upper()
    assert lower < 0 < upper < 1
    cost_floor = (1-upper)/(upper-lower)
    assert arb("7274.9177") < cost_floor < arb("7274.9179")
    print("公共牛顿矩与相等真空矩所需的负权总量下界", cost_floor.str(15))
    print("ΔΣgz、抵消系数、H 的全区间上下界",
          *[z.str(14) for z in (delta_moment, factor, lower, upper)])


if __name__ == "__main__":
    verify()
