"""A gauge-compatible competitor rules out cross-zero gain minimization as a CP selector."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
model = run_path(str(root/"ch39-four-contact-channels.py"))
gain_model = run_path(str(root/"ch39-shell-ckm-schur-gain.py"))
gauge = run_path(str(root/"ch39-schur-winner-gauge-yukawa.py"))
word = ("1242242122435512424315003425355335515353305521402552124342222044"
        "1115131112531513425051220343514220524242513002535021411502400211")


def verify():
    response = model["fixed"]["native"]["plane_response"](acb.zeta_zero(11618).imag)
    params = model["center_fit"]((".38251745472","3.40252006447","2.53367410385"),response)
    delta,jcp = model["phase_readout"](params,response)
    bound = arb(".38252")
    assert 0 < params[0] < bound and params[1] > 0
    assert arb(".16134358") < delta < arb(".16134360")
    print("11618 three-angle gain, phase and J",params[0].str(14),delta.str(14),jcp.str(14))
    cap = arb(".054")
    for z in acb.zeta_zeros(10815,1278):
        assert abs(gain_model["shell_coefficient"](z.imag)) < cap
    c2 = model["fixed"]["shell"]["axis_cos_sq"]
    factor = 2*(c2*(1-c2)).sqrt()*5/arb(867).sqrt()
    s,t,u = (model[k] for k in ("s","t","u"))
    a,b = (1-s*s).sqrt(),(1-t*t).sqrt()
    vtd = u + factor*bound*cap
    lower_cos = ((s*t)**2 + (a*b*u)**2-vtd**2)/(2*s*t*a*b*u)
    assert lower_cos.acos() < arb(".676")
    print("any candidate with gain <= .38252: absolute CP phase <",lower_cos.acos().upper().str(12))
    s,t,u = [v+arb(0,model["shell"]["sine_radius"].upper()) for v in (s,t,u)]
    a,b = (1-s*s).sqrt(),(1-t*t).sqrt()
    vtd = u+factor*bound*cap
    full_cos = ((s*t)**2+(a*b*u)**2-vtd**2)/(2*s*t*a*b*u)
    assert full_cos.lower().acos() < arb(".728")
    errors = gauge["gauge_match"](word,11618,"26.802773689811964")
    assert all(e < arb("1.1e-7") for e in errors)


if __name__ == "__main__":
    ctx.prec = 192
    verify()
