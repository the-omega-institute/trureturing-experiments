"""An explicit norm-five word preparing gauge coordinates for zero 11737."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

from flint import arb,ctx

root = Path(__file__).parent
words = run_path(str(root/"ch39-hurwitz-gauge-words.py"))
reference = run_path(str(root/"ch39-grid-gauge-readout.py"))
model, q = reference["model"],reference["q"]
WORD = ("0003403050530214353040414300511533135111112244303035533353331352"
        "2505143150415113552242135124244224124030334213341343040334153553")


def verify():
    ctx.prec = 144
    word = list(map(int,WORD))
    assert len(word) == 128 and all((a^1)!=b for a,b in zip(word,word[1:]))
    state = (1,0,0,0)
    for j in word:
        state = words["multiply"](state,words["GENERATORS"][j])
    assert sum(x*x for x in state) == 5**128
    a,*imag = state
    couplings = [abs(Q(v,a)) for v in imag]
    assert all(lo < g < hi for g,(lo,hi) in zip(couplings,words["BANDS"]))
    inputs = [q(g*g) for g in couplings] + [y*y for y in model["yukawa_centers"]]
    time = q(23) + q(10*(2*13349175109351+1))/2**46
    for _ in range(256):
        inputs = reference["advance"](inputs,time/256,14)
    values = [1/(u*r) for u,r in zip(inputs,model["ratios"])]
    target = model["target"](11737,time)
    residual = max(abs(x-target).upper() for x in values)
    assert residual < arb(2)**-23
    # The same full-box derivative bounds are replayed in the prior-preparation certificate.
    tolerance = arb(2)**-21
    max_shift = (2*tolerance+abs(values[2]-values[0]))/q(".12")
    time_band = time + arb(0,max_shift)
    color_band = values[2] + arb(0,tolerance+q(".16")*max_shift)
    assert model["target"](11736,time_band)>color_band
    assert model["target"](11738,time_band)<color_band
    print("128-step word, maximum normalized residual at 11737",residual.str(10),
          "full-window possible matching time displacement <",max_shift.upper().str(10))


if __name__ == "__main__":
    verify()
