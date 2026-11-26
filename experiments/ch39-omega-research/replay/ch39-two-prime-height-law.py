"""Drive both prime clocks with the same shell turn using the dyadic heights alone.

Prime rates are computed from finite rational sums at 1/3 and 1/2.
The phase columns from the original archive are read only by the checker.
"""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path
import json


def nearest(x):
    return (x + Q(1,2)) // 1


def finite_rate(prime, ticks=5):
    x = Q(prime-1,prime+1)
    term = x
    total = Q(0)
    for k in range(112):
        total += term/(2*k+1)
        term *= x*x
    return 2*ticks*total


def rate_numerators(ticks=5):
    root = Path(__file__).resolve().parent
    word = run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    # This native positive microstep is produced by the rational arc subtraction.
    exponent = -43235114615559591407865398948523
    (lo,hi), (low,high) = word["power"](exponent)
    unit = word["DEN"]
    assert lo*2 > unit and 0 < low <= high and high < unit >> 107
    middle = Q(low+high,2*unit)
    n,m = 1,0
    for _ in range(40):
        n,m = 4*n+17*m,n+4*m
    assert n*n-17*m*m==1
    coeff = [Q(nearest(finite_rate(p,ticks)*Q(m,n)/middle*2**24),2**24) for p in (2,3)]
    return [int(rate*2**24) for rate in coeff]


def main():
    root = Path(__file__).resolve().parent
    word = run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    record = json.loads((root/"ch39-five-tick-record.json").read_text())
    den = 1<<record["bits"]
    exponent = -43235114615559591407865398948523
    coeff = [Q(n, 2**24) for n in rate_numerators()]
    rows=[]
    for idx,height,_,_,_ in record["rows"]:
        # No input from phase columns or weights is used to construct these words.
        ticks=[exponent*nearest(rate*Q(height,den)) for rate in coeff]
        rows.append([idx,[0,ticks[0]],[0,ticks[1]]])
    word["check_output"](rows,record)
    print("Common generator C only, microstep exponent",exponent)
    print("Rate numerators over 2^24:",[int(rate*2**24) for rate in coeff])
    return rows


if __name__ == "__main__":
    main()

