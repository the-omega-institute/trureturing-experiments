"""Calibrate the additional 292 heights, adjacent target, and exact positive weights."""
import json
from pathlib import Path
from runpy import run_path
import sys
from flint import acb, arb, arb_series, ctx

root = Path(__file__).resolve().parent


def additions(case):
    cut = case["cut"]
    ctx.prec,ctx.cap = 640,3
    floor = run_path(str(root/"ch39-tick-separators-origin.py"))["floor_at"]
    heights = [v.imag for v in acb.zeta_zeros(3201,cut-3200)]
    logs = [arb(z).log() for z in (2,3,arb(3)/2,6)]
    extra = [[floor(h,400),*[floor((h*z).cos(),72) for z in logs]] for h in heights]
    def xi(j):
        t = arb_series([arb(1)/2+j,1])
        return t*(t-1)/2*((-t/2)*arb.pi().log()).exp()*(t/2).gamma()*t.zeta()
    origin, high = xi(0),xi(1683)
    target = floor(high[1]/(2*1683*high[0]*(origin[2]/origin[0])),288)
    a,w = run_path(str(root/"ch39-extra-offset-cutoffs-origin.py"))["compute"](case)
    assert a == case["a"]
    return {"case":{**case,"weights":w},"tail":extra,"target1683":target}


if __name__ == "__main__":
    path = root/"ch39-adjacent-offset-certificate.json"
    if len(sys.argv)==2:
        record = additions(json.loads(Path(sys.argv[1]).read_text()))
        path.write_text(json.dumps(record)+"\n")
    else:
        record = json.loads(path.read_text())
        assert additions(record["case"]) == record
    print("Adjacent RH target, 292 height/phase rows, and strictly positive exact arithmetic fit certified.")
