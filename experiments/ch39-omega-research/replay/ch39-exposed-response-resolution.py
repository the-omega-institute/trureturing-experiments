"""用两条球壳弧、真实零点计数和规范均值给出有限分辨率。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

certificate = run_path(str(Path(__file__).with_name("ch39-exposed-flavor-response.py")))
study = certificate["study"]
geo = study["geo"]
center_index = 11431
slope = arb("1.59828332e-8")
cutoff = 1_000_000


def zero_count(height):
    count = arb(height).zeta_nzeros().unique_fmpz()
    assert count is not None
    return int(count)


def spectrum_gap(A, B, slope, curvature, allowance, *, individual=(), readout=None, response_at=None):
    height = acb.zeta_zero(center_index).imag
    offset = study["offset"]
    numerator = height**2 + offset
    def ratio(h):
        return numerator / (h*h + offset)

    a, b = geo["a"], geo["b"]
    pi = arb.pi()
    step = 2*pi/a
    if response_at is None:
        response_at = study["native"]["plane_response"]
    reference = response_at(height)
    def score(h):
        diff = response_at(h)-reference
        if readout is not None:
            return readout(diff)
        return (A*(diff[0,0]-diff[1,1]).real/2
                + B*diff[0,1].real - diff[0,1].imag)

    def check(index, h):
        if index == center_index or index in individual:
            return
        target = slope * ratio(h)
        angle = h-height
        lower = curvature * (4-2*(a*angle).cos()-2*(b*angle).cos())
        if not lower > target:
            assert score(h)+allowance > target, index

    start = 200
    for index, zero in enumerate(acb.zeta_zeros(1, zero_count(start)), 1):
        check(index, zero.imag)

    # Nearest first-phase return. Include an extra return beyond either end.
    first = int(((start-height)/step).floor().unique_fmpz())-1
    last = int(((cutoff-height)/step).ceil().unique_fmpz())+1
    assert height+(first+arb(1)/2)*step < start
    assert height+(last-arb(1)/2)*step > cutoff
    windows = 0
    candidates = 0
    for k in range(first, last+1):
        mid = height+k*step
        cell_min = max(arb(start), (mid-step/2).lower())
        target = slope*ratio(cell_min)/curvature
        # For |t| <= pi, 2-2*cos(t) >= (2*t/pi)^2.
        half_width = (pi/(2*a)*target.sqrt()).upper()
        window = arb(mid, half_width)
        if 2-2*(b*(window-height)).cos() > target:
            continue
        # Enlarge on both sides so counts also cover possible endpoint zeros.
        lo = max(start-1, int(window.lower().floor().unique_fmpz())-1)
        hi = min(cutoff+1, int(window.upper().ceil().unique_fmpz())+1)
        if hi <= lo:
            continue
        windows += 1
        left, right = zero_count(lo), zero_count(hi)
        candidates += right-left
        if right > left:
            for index, zero in enumerate(acb.zeta_zeros(left+1, right-left), left+1):
                check(index, zero.imag)
    print("Phase cells, remaining windows, zero indices in enlarged windows",
          last-first+1, windows, candidates)

    gap = slope*(1-ratio(arb(cutoff)))-allowance
    h2 = acb.zeta_zero(2).imag
    upper = score(h2)/ratio(h2)
    return gap,upper


def verify():
    ctx.prec = 192
    A,B,_ = certificate["polynomial"]()
    gap,upper = spectrum_gap(A,B,slope,certificate["curvature"],arb(0))
    assert gap > arb("1.59808e-8") and upper < arb("1.598284e-8")
    print("Exact gauge mean, no center: certified lower and asymptotic upper",
          gap.str(17), upper.str(17))

    # Opposite gains cancel even orders of the ordered evolution.
    gain = arb(2)**-20
    norm = arb(9)/4
    assert norm > 2/(arb("1.3")-4*arb.pi()/32)
    functional_norm = (1+A*A+B*B).sqrt()
    error = 2*functional_norm*norm*((gain*norm).cosh()-1)
    finite_gap = gap-error
    assert finite_gap > arb("1.59704e-8")
    print("Gain, finite-gain correction, finite-gain lower", gain.str(12),
          error.str(12), finite_gap.str(15))


if __name__ == "__main__":
    verify()
