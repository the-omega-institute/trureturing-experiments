"""共同非负固有时轮廓下，三个有限类空动量的超荷逆耦合组合。"""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    gravity = run_path(str(root / "ch39-chiral-induced-gravity-vacuum.py"))
    spectrum = gravity["spectrum"]
    model, _, _, _, _ = spectrum["torus"]["local"]["inputs"]()
    ctx.prec = 224
    bands = spectrum["shell_bands"](model)
    assert gravity["prefix_cap"](bands, 32) == 0
    # 超荷在规范矢量上为零。每个重费米层使用 Dirac 权。
    fermions = (3 + 2*16)*sum(weyl[0] for weyl in model["weyl"])
    scalars = (1+16)*model["higgs"][0]
    assert (3*sum(weyl[0] for weyl in model["weyl"])+model["higgs"][0]-model["beta"][0]).contains(0)
    pi = arb.pi()
    energy = arb(1000)
    p2 = (energy/arb("1.22089e19"))**2
    moment_bound = 2*pi*arb("1e-122")/31
    one_point_error = (48*fermions + 24*scalars)*moment_bound/(4*pi*p2**2)
    # 对 Y(4P)-5Y(2P)/4+Y(P)/4 三个余项逐项取界。
    three_point_error = one_point_error*(arb(1)/256 + arb(5)/64 + arb(1)/4)
    logarithmic = arb(3)/4*model["beta"][0]*arb(2).log()/(2*pi)
    inverse_error = (logarithmic-three_point_error)*arb(2)/5
    step_error = (logarithmic-three_point_error)*arb(4)/5
    assert three_point_error < arb("2.714e-57")
    assert inverse_error > arb(".1356") and step_error > arb(".2713")
    print("三个动量下共同轮廓的组合余差、标准轻场的同一组合绝对值",
          three_point_error.upper().str(14), logarithmic.str(14))
    print("至少一项逆超荷耦合或一段变化的必要绝对误差",
          inverse_error.lower().str(14), step_error.lower().str(14))

    negative_vacuum = (logarithmic/three_point_error-1)*arb("1e-122")/2
    assert negative_vacuum > arb("6.24e-67")
    print("允许翻号后，精确恢复三点对数变化所需的负真空贡献下界",
          negative_vacuum.lower().str(14))

    # 在普朗克质量以下还可容纳 51 次的局部接触多项式检验。
    degree = 51
    assert energy*2**(degree+2) < arb("1.22089e19")
    attenuation = Fraction(1)
    for j in range(1, degree+1):
        attenuation *= Fraction(4**j - 1, 4**j + 1)
    assert attenuation > Fraction(1, 2)
    log_step = model["beta"][0]*arb(2).log()/(2*pi)
    floor = log_step*arb(3)/20 - one_point_error
    assert floor > arb(".0678")
    print("任意次数多项式经相应取样后的最坏逆耦合差下界", floor.lower().str(15))

    # 在原模型的电弱对称点，轻 Higgs 完全无质量，I(0)=∫H。
    fields = run_path(str(root / "ch39-chiral-minimal-neutral-gravity.py"))["fields"]
    assert all(h > 0 for values in bands.values() for h in values)
    rows = fields({q: [h.lower() for h in values] for q, values in bands.items()})
    heavy_k_error = sum(abs(a)/z.lower() for _, a, z in rows[1:])*moment_bound/(2*pi*arb(1).exp())
    area_lower = 12*pi*(1-heavy_k_error)
    # 所有带超荷标量贡献的首项同号，只保留原生轻双重态。
    drop = 6*model["higgs"][0]*area_lower*arb(3)/(4*4*pi*p2)
    drop -= one_point_error*(1+arb(1)/16)
    assert drop > arb("1e32")
    print("共同轮廓从 1 TeV 到 2 TeV 的逆超荷耦合下降量下界",drop.lower().str(15))


if __name__ == "__main__":
    verify()
