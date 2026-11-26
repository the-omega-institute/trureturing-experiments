"""Uniform width bound on the spectral factors served by one fixed low-energy source."""
from pathlib import Path
from runpy import run_path

from flint import arb

flow = run_path(str(Path(__file__).with_name("ch39-ckm-gauge-displacement.py")))
gauge, coefficients, loop = flow["model"],flow["volume"]["running"],flow["loop"]


def verify():
    low, high, logs = flow["volume"]["input_bounds"]()
    def check(step, g, traces):
        rates = [(-2/(r*loop))*(beta+sum(b*x-d*y for b,x,d,y in
                 zip(row,g,yukawa,traces))/loop)
                 for r,beta,row,yukawa in zip(gauge["ratios"],coefficients["BETA"],
                                                 coefficients["B"],coefficients["D"])]
        assert rates[2]-rates[0] > arb(".13") and abs(rates[0]) < arb(".06")
        assert 1/(g[0]*gauge["ratios"][0]) > arb("2.5")
    gauge["evolve"]([a.union(b) for a,b in zip(low,high)], logs, arb(33), 8192,
                    diagonal_quarks=True,observe=check)
    drift = [a+b for a,b in flow["displacement"](33)]
    tolerance = arb(2)**-21
    width = 2*(2*tolerance+drift[0]+drift[2])/arb(".13")
    change = 2*(tolerance+drift[0])+arb(".06")*width
    floor = arb("2.5")-tolerance-drift[0]
    log_width = 2*width + change/floor
    assert log_width < arb(".000037")
    assert (log_width.exp()-1) < arb(".000037")
    assert (1-arb("1.5e-6"))/(1-arb(".000085")) > log_width.exp()
    print("Any two matching clocks differ by less than",width.upper().str(8))
    print("Maximum spectral-factor ratio below 1 +",log_width.expm1().upper().str(8))


if __name__ == "__main__":
    verify()
