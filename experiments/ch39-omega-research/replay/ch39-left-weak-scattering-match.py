"""在共同左弱交换及独立左手接触项下检验电子与 μ 子的绝对率反演。"""
from pathlib import Path
from runpy import run_path

from flint import arb

study = run_path(str(Path(__file__).with_name("ch39-absolute-electron-muon-rates.py")))
poles, light, q, window = (study[key] for key in ("poles", "light", "q", "window"))
channels, physical, anchor = (study[key] for key in ("channels", "physical", "anchor"))


def weak_readout(time, logs, targets):
    (exact_length, exact_z, rho, _, _, _), _, _ = poles["retune"]["predictions"](logs, time)
    m1, m2 = poles["reference_readings"](logs)[:2]
    mass_window = window(q("1e-10"))
    length, z, _ = poles["charged_geometry"](m1*mass_window, m2*mass_window)
    reference = physical*light["graphs"](logs, time)[0]*physical.transpose()
    a, b, d = [reference[i, j] for i, j in channels]
    assert min(a, b, d) > 0
    common_rate = window(q("2e-5"))
    charged_cap = q("1/4")
    raw_error = q("1e-9")
    root_error = raw_error/(1+(1-raw_error).sqrt())
    # 参考三重态交换在 [0,1/4] 内任取。独立局域匹配残差以绝对振幅计。
    mismatch = q("1e-9")
    assert a*((1-raw_error)*(1-q("2e-5"))).sqrt() > mismatch
    radius = root_error*(a+2*charged_cap)+mismatch/(1-q("2e-5")).sqrt()
    extract = study["extract"]
    weights = [abs(extract(length, z, *row)) for row in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    uncertainty = weights[0]*radius + root_error*(weights[1]*b+weights[2]*d)
    for sign in (1, -1):
        k = (extract(length, z, a, sign*b, d)+arb(0, uncertainty.upper()))*common_rate.sqrt()
        assert k > 0
        inferred = 15*length*(1/k-poles["offset"])
        print("weak readout", sign, inferred.lower().str(12), inferred.upper().str(12), flush=True)
        if sign == 1:
            assert targets[2] < inferred < targets[0] and inferred.contains(targets[1])
        else:
            assert inferred < 0

    sensitivity = study["extract"](exact_length, exact_z, q(1), q(0), q(0))
    assert sensitivity > 0 and study["extract"](exact_length, exact_z, a, b, d).overlaps(1/rho)
    neighbours = (targets[0], targets[2])
    witnesses = study["pair_witnesses"](logs, time, reference, neighbours)
    relative_windows = (("28.94", "29.03"), ("-15.64", "-15.58"))
    for target, candidate_ll, bounds in zip(neighbours, witnesses, relative_windows, strict=True):
        k = 1/(poles["offset"]+target/(15*exact_length))
        contact = (1/rho-k)/sensitivity
        assert (candidate_ll+contact).overlaps(a)
        ppb = contact/a*10**9
        assert q(bounds[0]) < ppb < q(bounds[1])
        print("LL matching shift", contact.str(12), "relative ppb", ppb.str(10), flush=True)


def verify():
    time = poles["study"]["nominal"][2]
    targets = [poles["study"]["comp"]["model"]["target"](n, time) for n in (11430, 11431, 11432)]
    for logs in poles["retune"]["log_sets"]:
        weak_readout(time, logs, targets)


if __name__ == "__main__":
    verify()
