"""A rational finite-time comparison independent of all eight initial Yukawas."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path


def bounds(ratio_max, w_ratio, a, b, c, time_cap):
    pi_lo, pi_hi = F(333,106), F(355,113)
    ulo,uhi = 1/(8*pi_hi**2),1/(8*pi_lo**2)
    zmin = F(739,1000)**2*(1+3/ratio_max**2)
    A = F(8192,15)-15*a-F(12,5)*b-F(54,5)*c
    assert A>0 and ratio_max**2<=w_ratio
    assert F(22169,30)*w_ratio>144
    rates = [ulo**2*A*zmin/2-2*uhi/3,
             ulo*(F(23,6)*a-F(54,5)-F(297,20)*w_ratio),
             ulo*(F(16,5)*b-F(72,5)),
             ulo*(F(44,5)*c-F(108,5))]
    speed = min(rates)
    assert speed>0 and 1/(speed*zmin)<time_cap
    return speed,1/(speed*zmin)


def verify():
    old = run_path(str(Path(__file__).with_name("ch39-gauge-yukawa-fixed-points.py")))
    old["model"]()
    for case in ((F(167,100),F(14,5),15,5,F(27,10),52),
                 (F(8,5),F(64,25),F(66,5),5,F(53,20),42)):
        speed, duration = bounds(*case)
        print("Scalar comparison coefficient >",round(float(speed)-1e-7,7),
              "; finite logarithmic duration <",case[-1])
        assert duration>0


if __name__ == "__main__":
    verify()
