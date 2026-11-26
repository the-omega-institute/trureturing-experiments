"""Bounds on the single held-out reading for any five-tick measure in the first 2323 zeros."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
base = run_path(str(root/"ch39-cutoff-thresholds.py"))
data,B = base["data"],base["Bounds"]
cheb = run_path(str(root/"ch39-tick-separators-integers.py"))["chebyshev_upper"]
OFFSETS = base["OFFSETS"]
HEIGHT = base["HEIGHT"]
GRID = 1 << 336


def time_term(j,h,c):
    # Upper bound, in units of 1/GRID, for c times the time kernel.
    z = (h+int(c>0))**2
    return -((-c*z*GRID)//(z+(j*HEIGHT)**2))


def main():
    lows = data("ch39-minimum-tick-finite-inputs.json")["levels"]
    old = [*lows[:43],data("ch39-fixed-register-extra-target.json")["floor"]]
    prefix = data("ch39-tick-separators.json")["prefix"][:2323]
    assert len(prefix)==2323 and prefix[-1][0]+1 < 2**415
    floor = data("ch39-adjacent-offset-certificate.json")["target1683"]
    new = B(Q(floor,2**288),Q(floor+1,2**288))
    ranges = []
    for witness in data("ch39-adjacent-range-data.json"):
        s = witness["direction"]
        unit = 1 << witness["depth"]
        coeff = witness["a"]
        phase = witness["c"]
        pd = witness["phase_den"]
        assert len(coeff)==len(old)==len(OFFSETS) and len(phase)==4
        assert s in (1,-1)
        max_error = None
        for h,*angles in prefix:
            e = time_term(1683,h,s*unit)
            e += sum(time_term(j,h,-a) for j,a in zip(OFFSETS,coeff))
            for a,angle in zip(phase,angles):
                l,r = cheb(angle,5)
                e -= (a*(l if a>=0 else r)*unit*GRID)//(pd*2**72)
            max_error = e if max_error is None else max(e,max_error)
        # The original word differs from the calibrated cosine by < 2^-96.
        slack = Q(max_error,unit*GRID)+Q(sum(map(abs,phase)),pd*2**96)
        target = sum(Q(v*(y+int(v>0 and j!=0)),unit*2**288) for v,y,j in zip(coeff,old,OFFSETS))
        allowance = (Q(sum(map(abs,coeff)),unit)+Q(sum(map(abs,phase)),pd))/2**272
        bound = target+slack+allowance
        weights = base["exact_fit"](witness,Q(1,100000))
        value = B(0)
        for i,w in weights.items():
            h = prefix[i-1][0]
            k = (1683*HEIGHT)**2
            value += w*B(Q(h*h,h*h+k),Q((h+1)**2,(h+1)**2+k))
        actual = value if s==1 else B(-value.hi,-value.lo)
        assert bound>=actual.lo and bound-actual.lo<Q(1,10**13)
        ranges.append(B(actual.lo,bound) if s==1 else B(-bound,-actual.lo))
    hi,lo = ranges
    assert Q(".0359539100577")<hi.lo<=hi.hi<Q(".0359539100579")
    assert Q(".0359538933428")<lo.lo<=lo.hi<Q(".0359538933431")
    assert new.lo-hi.hi>Q(1,2**20) and new.hi-lo.lo<Q(1,2**19)
    print("Attainable maxima/minima bracketed for the original 2323-height cutoff.")
    print("With the previous data exact or at 272 bits, 20 bits on the extra reading reject every realization; 19 accept all.")


if __name__=="__main__":
    run_path(str(root/"ch39-prime-clock-height-horizon.py"))["main"](5,15)
    main()
