"""Dyadic inverse-temperature thresholds for the two-port fermion addition response."""
from pathlib import Path
from runpy import run_path
from flint import arb, acb, ctx

finite = run_path(str(Path(__file__).with_name("ch39-finite-arithmetic-height.py")))
GRID = 2**64


def free_energy(t,beta):
    return t+((1+(-2*beta*t).exp()).log()-arb(2).log())/beta


def load(h,coherent):
    return h/(13+12*h) if coherent else h/208


def verify():
    ctx.prec = 256
    cases = [(10875,195923358475422178854846,[(False,60),(True,56)],5,arb(".015"),-1),
             (10854,195587219045243829165080,[(True,58)],7,arb(".06"),1)]
    for index,target,readouts,slope_floor,weight_slope,direction in cases:
        zero = acb.zeta_zero(index).imag
        best = finite["weight"](zero)
        def datum(n):
            g = arb(n)/GRID
            return finite["weight"](g),abs(finite["finite_value"](g))
        base_weight,base_t = datum(target)
        other_weight,other_t = datum(target+direction)
        assert other_t > base_t
        for coherent,threshold in readouts:
            f = 1/arb(13 if coherent else 208)
            small = 2**(threshold-1)
            diff = free_energy(other_t,small)-free_energy(base_t,small)
            assert other_weight/base_weight-1 > f*diff
            # The certified zero-temperature score outside the local cell has relative gap 1e-10.
            wide_gain = arb("1.02")  # covers W/W* on both outer bands
            envelope = finite["source"]["modulus"](-arb(59)/108)
            low = 10567 if index==10854 else 10603
            assert envelope/(arb(low)**2+arb(1)/4)/best<wide_gain
            level = base_weight/(1+load(base_t,coherent))/best
            soft = arb(2).log()/2**threshold
            assert 1-arb("1e-10")+wide_gain*f*soft < level
            coeff = 1/(13+12*arb(".02")) if coherent else 1/arb(208)
            radius = (1-level+level*coeff*(arb("3.65e-22")+soft))/(slope_floor*level*coeff-weight_slope)
            count = int((radius*GRID).ceil().unique_fmpz())+2
            assert count<64
            rows = {n:datum(n) for n in range(target-count,target+count+1)}
            for b in range(threshold,72):
                beta = 2**b
                winner = base_weight/(1+load(free_energy(base_t,beta),coherent))
                assert all(winner>w/(1+load(free_energy(t,beta),coherent))
                           for n,(w,t) in rows.items() if n!=target)
            # Enclose all beta >= 2^72 at once, without testing temperatures individually.
            spread = arb(2).log()/2**72
            floor = base_weight/(1+load(base_t,coherent))
            assert all(floor>w/(1+load(max(arb(0),t-spread),coherent))
                       for n,(w,t) in rows.items() if n!=target)
            print(index, "coherent" if coherent else "trace",
                  "first dyadic inverse temperature exponent",threshold,"local grid radius",count)


if __name__ == "__main__":
    verify()
