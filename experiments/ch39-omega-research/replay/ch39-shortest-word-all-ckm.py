"""Certify the earliest true-zero match for arbitrary unitary quark mixing."""
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

root = Path(__file__).parent
s = run_path(str(root/"ch39-shortest-matching-rotor.py"))
mixing = run_path(str(root/"ch39-ckm-gauge-displacement.py"))
q = s["q"]


def verify():
    ctx.prec = 128
    # The broad-window exclusion and 1278-label envelope already allow any
    # opposite-family overlap between zero and the entire opposite trace.
    enclosure = s["Enclosure"]()
    assert s["anchor"] + s["volume"]["radius"] < 27
    change = [a+b for a,b in mixing["displacement"](27)]
    assert all(v < q(limit) for v,limit in zip(change,("6.24e-8","7.034e-8","1.002e-7")))
    enclosure.center_error = [e+arb(0,v) for e,v in zip(enclosure.center_error,change)]
    s["exclude_other_orbits"](enclosure)
    a,*imag = s["match_coord"]
    ys = s["model"]["yukawa_centers"].copy()
    ys[2] *= q(".9936504245091254")
    time = q("26.8162148027")
    state = [(arb(i)/a)**2 for i in imag]+[v*v for v in ys]
    for _ in range(128):
        state = s["reference"]["advance"](state,time/128,14)
    target = s["model"]["target"](s["match_label"],time)
    residuals = [abs(1/(r*u)-target)+d for r,u,d in zip(s["ratios"],state,change)]
    assert all(v < q("4.693e-7") for v in residuals)
    print("All unitary mixing paths, unique shortest label 11439, residual bounds",
          *[v.upper().str(9) for v in residuals])


if __name__ == "__main__":
    verify()
