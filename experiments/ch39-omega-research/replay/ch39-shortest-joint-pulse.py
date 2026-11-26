"""Compare one through eight clock steps with the same local channel graph and inputs."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

local = run_path(str(Path(__file__).with_name("ch39-local-product-contact.py")))
contact = local["contact"]


def verify():
    ctx.prec = 192
    spectrum = local["graph_spectrum"]()
    windows = contact["clock_windows"](range(1,9))
    centers = contact["flavor"]["shell"]["centers"]
    tolerances = [arb(".000055")]*3+[arb(".006")]
    for length,rows in windows.items():
        n,energy,effective,excitation,probability,drift,vals,residual = contact["pulse_measurement"](
            rows,spectrum,length)
        scaled = [abs(v-c)/t for v,c,t in zip(vals,centers,tolerances)]
        if length <= 5:
            assert scaled[0] > 1
        else:
            assert n==11431 and all(e<1 for e in scaled) and residual < arb(2)**-21
        if length==4:
            assert n==11431 and scaled[0]>arb("2.089") and scaled[1]>arb("1.115")
        if length==6:
            assert probability>arb(".999832") and residual<arb("3.97e-7")
        print(length,n,probability.str(8), "flavor errors",*(v.str(6) for v in scaled),
              "gauge drift",drift.str(9),"sufficient residual",residual.str(7))


if __name__=="__main__":
    verify()
