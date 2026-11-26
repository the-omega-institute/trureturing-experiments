"""Compare finite-volume scores with a measured Behrend-shell four-tap readout."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    ctx.prec = 112
    degree = 32
    cutoff = arb(".070")
    a = (degree*(arb(".0795438")/cutoff).acosh()).cosh()
    n = (degree*(arb(".079546")/cutoff).acosh()).cosh()
    leakage = 1/(arb(".001")*a)**2
    floor = (a/n)**2*arb(".001")**2*(1-leakage)
    rounding = 2/(2**45*floor)
    retained = ((1-leakage)*(1-rounding)/((1+leakage)*(1+rounding)))**2
    per_coefficient = 4*(1-retained)*(arb(13)/4)**2
    root_two = arb(2).sqrt()
    filter_error = per_coefficient*3*(1+root_two)
    # The rational all-pass response preserves |D|=sqrt(2). On its closed
    # disk this bounds the derivative of the ideal polynomial (18,-19,-3,4)/13.
    derivative = (19+6*root_two+24)/13
    # Rounding the shell tangent to k/2^48 gives the input height 2^47/k.
    # The existing all-pass synthesis is certified over the entire rounded band.
    band_top = arb(11644)
    height_error = band_top**2/(arb(2)**48-band_top)
    assert acb.zeta_zero(10815).imag - height_error > 10570
    assert acb.zeta_zero(12092).imag + height_error < band_top
    phase_error = root_two*(arb(".00000052") + arb(2).log()*height_error)*derivative
    max_modulus = (18 + 19*root_two + 6 + 8*root_two)/13
    denominator_error = max_modulus*2*band_top*height_error/10000**2
    error = filter_error + phase_error + denominator_error
    assert 0 < error < arb(".000097")
    pointwise = run_path(str(Path(__file__).with_name("ch39-pointwise-volume-winner.py")))
    pointwise["verify"](error)
    print("64-degree shell readout, weight error times 10000^2 <",error.upper().str(10))


if __name__ == "__main__":
    verify()
