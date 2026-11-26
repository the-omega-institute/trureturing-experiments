"""A single independent RH time reading rejects the entire 50-site conditional family."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
s=run_path(str(root/"ch39-fixed-scale-common-reading.py"))
base,B,rounded=s["base"],s["B"],s["rounded"]


def extrapolation_cost(tri):
    rec=base["data"]("ch39-cutoff-threshold-inverses.json")[2]
    mat,_=base["observations"](rec["support"])
    inv=rec["inverse"]
    den=base["SCALE"]*2**rec["bits"]
    def apply(col,cells):
        lower=upper=0
        for v,(l,h) in zip(col,cells):
            lower+=v*(l if v>=0 else h)
            upper+=v*(h if v>=0 else l)
        return lower,upper
    error=0
    columns=list(zip(*inv))
    for i,row in enumerate(mat):
        size=0
        for j,col in enumerate(columns):
            lower,upper=apply(col,row)
            diag=den*int(i==j)
            size+=max(abs(lower-diag),abs(upper-diag))
        error=max(error,size)
    assert 2*error<den
    cells=[]
    for i in rec["support"]:
        r=tri.poles[i]
        v=r/(r+B(43**2))
        cells.append((v.lo*base["SCALE"]//1,-((-v.hi*base["SCALE"])//1)))
    size=sum(max(map(abs,apply(col,cells))) for col in columns)
    return Q(size,den-error)


def main():
    tri=s["ThreeFits"]()
    level=base["data"]("ch39-minimum-tick-finite-inputs.json")["levels"][43]
    target=B(Q(level,2**288),Q(level+1,2**288))
    prediction=sum((w*(tri.poles[i]/(tri.poles[i]+B(43**2))) for i,w in tri.fits[2].items()),B(0))
    error=prediction-target
    assert Q("4.1075e-50")<error.lo<=error.hi<Q("4.1077e-50")
    for k in range(2):
        shift=B(-43**2)*tri.difference(k,-Q(43**2))
        assert shift.lo>0
    noise=extrapolation_cost(tri)/2**272
    print("Original 48 weights, excess",float(error.lo),float(error.hi),"worst change from old tolerances",float(noise))
    assert error.lo-noise>Q(1,2**165) and error.hi<Q(1,2**164)
    print("First next ordinary offset j=43 rejects every 50-site model at 165 bits; the base fit passes at 164.")


if __name__=="__main__":
    main()
