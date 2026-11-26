"""Attain a competing continuous height using the diagonal-quark two-loop submodel."""
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx


def verify():
    ctx.prec = 144
    gauge = run_path(str(Path(__file__).with_name("ch39-floating-scale-winner.py")))
    q = gauge["q"]
    old_time = q("26.859794675")
    height = q("10617.645332791202380").union(q("10617.645332791202395"))
    index_height = acb.zeta_zero(10871).imag
    offset = q(gauge["fixed"]["OFFSET"])
    common = gauge["target"](10871,old_time)
    time = old_time + ((index_height**2+offset)/(height**2+offset)).log()/2
    assert 23 < time < 33
    logs = list(map(q, ("-25.4446915838 -12.9928374813 -1.5528331663 -23.8189018241 "
                        "-17.8423094899 -10.1208976879 -25.6505068564 -14.9499484326 "
                        "-9.2850641599").split()))
    start, yukawa, _ = gauge["evolve"]([common*r for r in gauge["ratios"]], logs,
                                        -time, 131072, diagonal_quarks=True)
    factors = [1/(x.sqrt()*c) for x,c in zip(start,gauge["centers"])]
    assert all(q(".999") < x < q("1.001") for x in factors)
    assert all(q(".99") < (y/2).exp()/c < q("1.01") for y,c in zip(yukawa,gauge["yukawa_centers"]))
    print("Feasible stationary competitor interval, weak and color input ratios",*[x.str(10) for x in factors[1:]])


if __name__ == "__main__":
    verify()
