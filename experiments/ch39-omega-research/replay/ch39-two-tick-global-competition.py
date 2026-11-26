"""Compare the radial two-slot clock with actual positive zeros outside the gauge window."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

root = Path(__file__).parent
pulses = run_path(str(root/"ch39-finite-port-pulses.py"))
independence = run_path(str(root/"ch39-zero-haar-clocks.py"))


def score(n):
    generator = pulses["generator"](acb.zeta_zero(n).imag)
    first, delayed = pulses["port_outputs"](generator,2,impulse=True)
    return abs((-7*first[0]+24*delayed[0])/25)**2


def verify():
    reference = score(11431)
    assert score(873)-reference > arb(".00000402")
    witness = score(2994)
    assert witness-reference > arb(".00000572")
    # Phase changes enter only the off-diagonal edge factors. Each 40x40
    # entry changes by at most 2*h when each of the two angles moves by h.
    radius = arb(2)**-24
    amplitude_error = (7*2+24*4)*80*radius/25
    assert witness.sqrt()-amplitude_error > 0
    lower_score = (witness.sqrt()-amplitude_error)**2
    assert lower_score > reference+arb(".000002")
    density = (radius/arb.pi())**2
    assert density > arb("3.59e-16")
    print("scores at 11431, 873, 2994",*(p.str(13) for p in (reference,score(873),witness)))
    print("certified torus square: excess above 2e-6, RH lower density above 3.59e-16.")


if __name__ == "__main__":
    verify()
