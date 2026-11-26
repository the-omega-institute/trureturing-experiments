"""Calibrate the single extra arithmetic target with an outward enclosure."""
from pathlib import Path
import json
from flint import arb, arb_series, ctx


def main():
    p=Path(__file__).with_name("ch39-fixed-register-extra-target.json")
    doc=json.loads(p.read_text())
    ctx.prec,ctx.cap=640,3
    def xi(j):
        s=arb_series([arb(1)/2+j,1])
        return s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
    base=xi(0); value=xi(doc["j"])
    reading=value[1]/(2*doc["j"]*value[0]*(base[2]/base[0]))
    assert int((reading*2**doc["bits"]).floor().unique_fmpz())==doc["floor"]
    print("Extra target verified:",reading.str(19))


if __name__=="__main__":
    main()

