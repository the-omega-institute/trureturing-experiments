"""由三条复谐振链的真实零点能检验端口计次选择。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

chain = run_path(str(Path(__file__).with_name("ch39-local-winding-chain.py")))
angle = chain["joint"]["winding"]["angle"]


def oscillator_energy(p):
    return 3*((arb.pi()/(4*p)).cot()-1)


def phase_energy(p, theta, radial):
    spacing = 2*arb.pi()/3
    quotient = p*theta/spacing
    base = int(quotient.floor().unique_fmpz())
    distance = min(abs(p*theta-j*spacing) for j in range(base-2,base+4) if j%3)
    value = 2*(1-distance.cos())/p
    return 3*value/(1+radial*value/2)


def verify():
    ctx.prec = 128
    target = 23226
    theta = angle(11431)
    radial = (arb(1)/6).union(arb(1)/4)
    lower_step, upper_step = 12/arb.pi(), 3*arb.pi()/2
    assert arb("3e-18")*lower_step>arb("1e-17")
    assert arb("4e-17")*upper_step<arb("2e-16")
    assert oscillator_energy(1).contains(0)
    current_phase = phase_energy(target,theta,radial)
    lower_candidate = 4*target+1
    min_weight = (current_phase-phase_energy(lower_candidate,theta,radial))/(oscillator_energy(lower_candidate)-oscillator_energy(target))
    max_weight = (phase_energy(target-1,theta,radial)-current_phase)/(oscillator_energy(target)-oscillator_energy(target-1))
    assert min_weight>0 and max_weight>min_weight
    assert min_weight > arb("3.0636e-19") and max_weight < arb("4.3433e-13")
    print("necessary normalized oscillator unit",min_weight.str(12),max_weight.str(12),flush=True)
    minimum_mode = 2*(arb.pi()/(2*target)).sin()
    assert max_weight*minimum_mode < arb("5.875e-17")
    print("necessary chain gap / phase energy unit, upper", (max_weight*minimum_mode).upper().str(12),flush=True)


if __name__=="__main__":
    verify()
