"""Integer certificate for both prime clocks on a full finite height range."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
from runpy import run_path


def main(ticks=5, height_bits=36):
    root = Path(__file__).resolve().parent
    words = run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    law = run_path(str(root/"ch39-two-prime-height-law.py"))
    unit = words["DEN"]
    (amin,_),(bmin,bmax) = words["power"](-43235114615559591407865398948523)
    assert amin*2 > unit and bmin > 0
    den = 1 << 320
    k = isqrt(17*den**2)
    assert k*k < 17*den**2 < (k+1)**2
    low = Q(k*bmin,den*unit)
    sine_hi = Q((k+1)*bmax,den*unit)
    assert sine_hi < Q(1,2)
    # On the positive branch, z <= arcsin(z) <= z+z^3 (0<=z<=1/2).
    # This follows from 1 <= (1-t^2)^(-1/2) <= 1+2*t^2 there.
    high = sine_hi + sine_hi**3
    errors=[]
    for p, num in zip((2,3), law["rate_numerators"](ticks)):
        rate = Q(num, 1<<24)
        x = Q(p-1,p+1)
        lower = law["finite_rate"](p,ticks)
        # Finite geometric remainder integrated from 0 to x for 10/(1-t^2).
        upper = lower + 2*ticks*x**225/(225*(1-x*x))
        errors.append(max(abs(rate*low-upper),abs(rate*high-lower)))
    # Two rounding errors of at most half a microstep each and the error
    # incurred by recording the height at 288 binary places.
    bound = high + 2**height_bits*sum(errors) + Q(4*ticks,2**288)
    assert bound < Q(773,1000*2**96)
    print("At tick",ticks,"and heights through 2^",height_bits,"all four sine/cosine errors are below 0.773 * 2^-96.")


if __name__ == "__main__":
    main()

