"""Explicit power guaranteeing a positive general quartic and tree J=0 margins."""
from pathlib import Path
from runpy import run_path
from flint import arb

base = run_path(str(Path(__file__).with_name("ch39-shell-gw-full-inert-potential.py")))
q=base["q"]


def check(m):
    gain=q(81,977)**m
    T,ZF,scalar,VM,raw,time=base["check"](gain)
    # Heavy doublet Yukawa indices alone obey d||Y_H||/dx <= 6 T ||Y_H||.
    heavy=raw.sqrt()*gain*(6*T*time).exp()
    Y=(12*T).sqrt()
    YH=arb(12).sqrt()*heavy
    # E is the Frobenius norm of the canonically normalized real quartic tensor difference.
    # The scalar reference sees the actual g_i and actual light Yukawa traces.
    hard=q(".025")
    L=12*scalar+3*hard+4*Y**2
    force=(8*Y*ZF+96*Y**3)*YH
    error=force*((L*time).exp()-1)/L
    print("power",m,"error",error.upper().str(8),"positivity budget",(6*VM).lower().str(8),
          "baseline Frobenius upper",ZF.upper().str(8),"time",time.upper().str(8))
    assert error<hard and error<6*VM
    assert scalar+error/2<8*arb.pi() and 72*T<8*arb.pi()
    return VM-error/6,scalar+error/2


if __name__=="__main__":
    margin,scatter=check(5)
    print("full-interval quartic lower coefficient above",margin.lower().str(8),
          "contact amplitude norm below",scatter.upper().str(8))
