"""Read the gauge target by comparisons with a normalized arithmetic reference."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, arb

leads = run_path(str(Path(__file__).with_name("ch39-unit-lead-gauge-shell.py")))


def verify():
    floating = leads["floating"]
    target = floating["gauge"]["shortest"]["model"]["target"]
    first_weight = 1/(acb.zeta_zero(1).imag**2+arb(1)/4)
    weight_sum = 1+arb.const_euler()/2-(4*arb.pi()).log()/2
    probability_floor = first_weight/weight_sum
    assert probability_floor > arb(".2164")
    c = arb(452)/15
    alpha = c*first_weight
    scale = 2**96
    code = (alpha*scale).floor().unique_fmpz()
    assert code is not None and 0 < code < scale
    coin = arb(code)/scale
    qmax = target(1,arb(23))
    lowest_rail = target(1,arb(33))*(1+alpha)/alpha
    assert lowest_rail > 57
    trial_scale = arb("26.81676714125598")
    assert target(11431,trial_scale)+arb("1.005e-7")+arb(".09")*(33-23) < 5
    rounding_error = qmax*(alpha-coin)/alpha
    truncation_error = qmax*(1-probability_floor/2)**512
    error = arb("1.3e-16")
    assert 0 <= rounding_error and rounding_error+truncation_error < error
    # A finite-only trial budget uses p_1 >= 1/N and (1-1/(2N))^(2N) < 1/2.
    assert qmax/arb(2)**85+rounding_error < error
    # Check the first-decision probability by dynamic programming on a finite
    # reference and compare with the diagonal compression formula exactly.
    weights = [Fraction(4,9),Fraction(1,3),Fraction(2,9)]
    toss = Fraction(int(code),scale)
    for p in weights:
        accept,calibrate = toss*p/2,weights[0]/2
        active,output = Fraction(1),Fraction(0)
        for _ in range(5):
            output += active*accept
            active *= 1-accept-calibrate
        assert output == accept/(accept+calibrate)*(1-(1-accept-calibrate)**5)
    # In a passive divider the candidate conductance falls with the zero weight.
    # Rayleigh monotonicity then favors 14840 over 11431 for all scales.
    free = floating["penalty"]["doublet"]["joint"]["contact"]["clock"]["free_ground"]
    score = floating["penalty"]["competition"]["score"]
    assert acb.zeta_zero(14840).imag > acb.zeta_zero(11431).imag
    assert free(score(11431),2)-free(score(14840),2) > arb(".00000874")
    leads["verify"](target_error=error)
    print("dyadic coin numerator",code,"reference first-label probability above",
          probability_floor.lower().str(9),"uniform arithmetic terminal error below",
          (rounding_error+truncation_error).upper().str(9))


if __name__ == "__main__":
    verify()
