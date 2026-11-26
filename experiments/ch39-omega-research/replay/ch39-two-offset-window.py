"""Uniform finite product probes between the squared offsets 1682^2 and 1683^2."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
model = run_path(str(root/"ch39-same-offset-square.py"))
base, B = model["base"], model["base"]["Bounds"]
A = 1682**2
D = 1683**2-A


def cell(n):
    return B(Q(n,2**288),Q(n+1,2**288))


def environment():
    data = base["data"]
    prefix = data("ch39-tick-separators.json")["prefix"]+data("ch39-adjacent-offset-certificate.json")["tail"]
    squared = data("ch39-same-offset-square-data.json")
    square = squared["square_floor"]
    first = cell(data("ch39-fixed-register-extra-target.json")["floor"])
    second = cell(data("ch39-adjacent-offset-certificate.json")["target1683"])
    targets = [cell(square),(B(A+D)*second-B(A)*first)/B(D)]
    return prefix,square,targets


def old_system(ids):
    prefix,square,_ = environment()
    rows,rhs = model["square_system"](ids,prefix,square)
    return rows[:44]+rows[45:],rhs[:44]+rhs[45:]


def main():
    data = base["data"]
    prefix,square,targets = environment()
    unit = 2**288
    err = Q(1,2**272)
    # Linear interpolation minus the true product is nonnegative and at most bend.
    # For each z, K_A^2 (D z/(1+A z))^2 <= (D/A)^2/16.
    bend = Q(D*D,64*A*A)
    assert bend<Q("0.000000023")
    lows = data("ch39-minimum-tick-finite-inputs.json")["levels"][:43]
    first = data("ch39-fixed-register-extra-target.json")["floor"]
    second = data("ch39-adjacent-offset-certificate.json")["target1683"]
    old = [*lows,first]
    sq_case = data("ch39-same-offset-square-data.json")["case"]
    adj_case = data("ch39-adjacent-offset-certificate.json")["case"]
    for case,values,square_case in ((sq_case,old+[square],True),(adj_case,old+[second],False)):
        coeff = [Q(c,2**40) for c in case["a"]]
        angles = [Q(c,2**case["depth"]) for c in case["c"]]
        assert coeff[-1]>0
        truth = sum(q*Q(y+int(q<0),unit) for q,y in zip(coeff,values))
        # Both verified separators are < 0.000001 on the first 3491 columns.
        old_error = err*(sum(map(abs,coeff[:-1]))+sum(map(abs,angles)))
        gap = (truth-Q(".000001")-old_error)/coeff[-1]
        if not square_case:
            gap = (gap*(A+D)-A*err)/D
        assert gap>bend+Q(".000005")
    for rec in data("ch39-two-offset-window-data.json"):
        ids = rec["support"]
        fit = base["exact_fit"](rec,Q(1,100000),old_system(ids),48,max(ids))
        assert max(ids)<=3492
        endpoints = [B(0),B(0)]
        for i,w in fit.items():
            h = prefix[i-1][0]
            z = B(Q(2**800,(h+1)**2),Q(2**800,h*h))
            k = B(1)/(B(1)+B(A)*z)
            for j,d in enumerate((0,D)):
                endpoints[j] += w*k/(B(1)+B(A+d)*z)
        direction = rec["direction"]
        margins = [v-t if direction==1 else t-v for v,t in zip(endpoints,targets)]
        print(direction,"endpoint margins",[float(m.lo) for m in margins])
        assert all(m.lo>bend for m in margins)
    assert sorted(rec["direction"] for rec in data("ch39-two-offset-window-data.json"))==[-1,1]
    # Bit threshold on a single, unscaled response once the old measurements are exact.
    first_cell,second_cell = cell(first),cell(second)
    hi,lo = Q(".0359539100579"),Q(".0359538933428")
    factor = Q(A+D,A)
    lower = first_cell.lo-hi-factor*(first_cell.hi-second_cell.lo)
    upper = factor*(first_cell.hi-lo)-(first_cell.lo-second_cell.hi)
    assert Q(1,2**20)<Q(".000001029")<lower<upper<Q(".000001090")<Q(1,2**19)
    print("The exact minimal highest index is 3492 on the full squared-offset window.")
    print("At step 3365/2^k the first resolving absolute-error exponent for cutoff 2323 is 20+k.")


if __name__=="__main__":
    main()
