"""A four-step shell operator acting uniformly on the two certified arithmetic spectra."""
from pathlib import Path
from runpy import run_path
from flint import arb

inert = run_path(str(Path(__file__).with_name("ch39-shell-gw-two-inert-heavy.py")))
flow = inert["flow"]


def T4(x):
    return 8*x**4-8*x**2+1


def check(index):
    heavy, data = inert["certify"](index)
    model = flow["full"]["model"]
    if index == 43:
        s,f,_,_ = model["m"]["transported_spectra"]()
    else:
        s,f,_ = model["scalar_models"]()
    gaps = flow["full"]["all_levels"](s,f)
    assert 0<gaps[0]<gaps[1]<1<gaps[2] and gaps[-1]<7
    stop = 1/T4(arb(4)/3)
    assert stop < arb(".083")
    actual = [abs(T4((4-x)/3))*stop for x in gaps[2:]]
    # Impossibility for even the discrete spectrum by interpolation at four
    # actual (simple or clustered) active-band nodes.
    nodes = [gaps[i] for i in (2,10,34,52)]
    norm = arb(0)
    for j,s in enumerate(nodes):
        lagrange = arb(1)
        for k,t in enumerate(nodes):
            if j!=k:
                lagrange *= t/(t-s)
        norm += abs(lagrange)
    lower = 1/norm
    assert lower > arb(".15")
    F,D = flow["fit"]["readouts"](index)
    total = sum(inert["hs"](M,M).real for M in (F,D))
    # All active heavy columns arrive no earlier than the third extra mass.
    Q,_ = flow["masses"](index)
    light = sum(inert["hs"](m,m).real for m in data)
    early = (heavy[2]/heavy[0]).log()/flow["loop"]
    assert 10*early*light<1
    cap = light/(1-10*early*light)+inert["g"]**2*total*max(actual)**2
    late = (Q/heavy[2]).log()/flow["loop"]
    assert 10*late*cap<1
    cap /= 1-10*late*cap
    assert cap < arb("2.6")
    print("zero",index,"lower on degree <=3 absolute gain",lower.lower().str(9),
          "max active gain for degree four",max(actual).upper().str(9),
          "total running Yukawa squared norm below",cap.upper().str(8))


if __name__=="__main__":
    for n in (43,195):
        check(n)
