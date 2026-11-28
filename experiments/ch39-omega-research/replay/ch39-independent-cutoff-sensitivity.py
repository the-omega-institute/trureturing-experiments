"""Independent logarithmic motion of signed cutoffs on the three visible spectra."""

from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def prefix_deficits(rows):
    """Upper and lower cumulative polarization, allowing all interval mass orderings."""
    ends = [end for _, _, z in rows for end in (z.lower(), z.upper())]
    lower, upper = [], []
    for t in ends:
        lower.append(-sum(g for g, _, z in rows
                          if (z.lower() if g > 0 else z.upper()) <= t))
        upper.append(-sum(g for g, _, z in rows
                          if (z.upper() if g > 0 else z.lower()) <= t))
    return min(lower), max(upper)


def check_spectrum(rows):
    light = [(g, a) for g, a, z in rows if z.is_zero()]
    heavy = [(g, a, z) for g, a, z in rows if not z.is_zero()]
    pi = arb.pi()
    c = arb("1e-122")
    assert len(light) == 1 and light[0][0] == -62
    assert light[0][1].contains(arb(1)/6)
    assert prefix_deficits(rows) == (62, 482)
    g_mass = sum(abs(g) for g, _, _ in heavy)
    a_mass = sum(abs(a) for _, a, _ in heavy)
    assert g_mass == 1544 and a_mass < 199
    assert sum(abs(a) for _, a, _ in rows) < 64*pi

    gap = min(z.lower() for _, _, z in heavy)
    charge = min(z.upper() for g, _, z in heavy if g == -24)
    vector = [z for g, _, z in heavy if g > 0]
    assert gap > arb("1.5e-11") and gap > 2000*c
    assert len(vector) == 18 and min(z.lower() for z in vector) > 2*charge
    scale = arb("2.435e18")/(724**2*arb("1.22089e19"))
    assert charge < (57*scale)**2 and min(z.lower() for z in vector) > (96*scale)**2

    # Enclose k/x, w/x², Dk/x, Dw/x² throughout x <= 16c.
    # The massive corrections are bounded by their coefficient mass times e^-125.
    u = arb(125)
    suppression = (-u).exp()
    k0 = 1/(12*pi)
    w0 = 62/(8*pi)
    ka = arb(0, a_mass*suppression/(2*pi))
    ga = arb(0, g_mass*suppression/(4*pi))
    k, dk = k0+ka, k0+ka
    w, dw = w0+ga, 2*w0+ga
    ratio_slope = (dw*k-w*dk)/(k*k)
    assert k > 0 and k < 2*k0 and 90 < ratio_slope < 96

    # P=(Dw-2cDk)/x, c/x <= u_j/2000 for every massive row.
    p_error = (g_mass*(1+u)+4*a_mass*u*u/2000)*suppression/(4*pi)
    p_slope = 62/(4*pi)+arb(0, p_error)
    assert 4 < p_slope < 6
    assert dw/186-2*dk < 0 and 2*dw/93-2*dk > 0

    # Lower bound each charged term in the remainder using E_2 <= E_1.
    for g, alpha, _ in heavy:
        if g < 0:
            assert g+4*arb(93)/2*alpha > 0
    assert all(alpha.contains(2) for g, alpha, _ in heavy if g == -24)
    assert all(g == 36 and alpha.contains(-6) for g, alpha, _ in heavy if g > 0)
    assert (4*arb(93)/2*2-24)/4 >= 87
    assert 18*(36+4*186*6) == 81000
    assert 87 > 81000*suppression
    a, b = 62/(4*pi), 958/(4*pi)
    assert (b/2+6)/(a-4) < 48 and 2*k0*96/4 < 50

    budget = arb(10)**6
    time = 2*c*(32**2)*budget/a
    drift = 2*time
    assert time < arb("4.1511e-114") and drift < arb("8.3022e-114")
    return gap, min(z.lower() for z in vector)/charge, time, drift


def main():
    model = run_path(str(Path(__file__).with_name("ch39-cutoff-scale-stationarity.py")))
    spectra = model["spectra"]()
    ctx.prec = 256
    for number, rows in spectra.items():
        gap, separation, time, drift = check_spectrum(rows)
        print(number, "deficit interval [62,482]; certified lower bounds", gap.lower(), separation.lower())
    print("T <= 10^6: maximum log(x) displacement <", time.upper().str(7),
          "Newton drift <", drift.upper().str(7))


if __name__ == "__main__":
    main()
