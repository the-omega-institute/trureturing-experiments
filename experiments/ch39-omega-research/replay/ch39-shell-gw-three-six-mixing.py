"""Test unrestricted relative mixing when one sector has exactly three modes."""
from itertools import product
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

seven = run_path(str(Path(__file__).with_name("ch39-shell-gw-seven-left-gram.py")))
model = seven["model"]
adj = model["adj"]
rows = seven["rows"]
count = seven["count"]


def orbit_frame(matrices):
    four = model["bands"]["four"]
    b0, b1 = map(model["array"], matrices[:2])
    left = four["unitary_trial"](np.linalg.eig(b1 @ np.linalg.inv(b0))[1])
    right = four["unitary_trial"](b0.conj().T @ left)
    return four["exact_frame"](left), four["exact_frame"](right)


def three_mode_bounds(matrices):
    diagonals = acb_mat([[b[i, i] for b in matrices] for i in range(3)])
    inv = diagonals.inv()
    for i, j in product(range(3), repeat=2):
        if i == j:
            continue
        entries = [b[i, j] for b in matrices]
        assert sum(abs(z)**2 for z in entries).sqrt() < arb(".00130")
        coeff = acb_mat([entries]) * inv
        assert abs(coeff[0, 3-i-j]) < arb(".000388")
        assert max(abs(z) for z in coeff.entries()) < arb(".00088")

    tol = arb("1e-6")
    rel = (1 + tol)/(1 - tol)
    for light, middle, heavy in [(".00000615", ".00311", ".8616"),
                                 (".0000135", ".000268", ".01401")]:
        small = arb(light) / arb(heavy) * rel
        low = arb(middle) / arb(heavy) / rel
        high = low * rel**2
        d = arb(".00131")
        assert arb(".00130")*(1 + 2*high**2).sqrt() < d
        top = arb(2).sqrt() * d / (1-high-4*d)
        other = top**2 + high
        cross = (arb(".000388") + 2*arb(".00088")*other + top**2 + small) / low
        # The largest off-pivot coordinate of the second left vector determines
        # its axis. Bound the opposite right coordinate, then the opposite left.
        large = ((1 - top**2)/2).sqrt()
        minor = cross / (1-top**2-(cross/large)**2).sqrt()
        second = (top**2 + minor**2).sqrt()
        first = (top**2 + second**2).sqrt()
        chord_sq = lambda t: 2*(1-(1-t**2).sqrt())
        frame_error = sum(chord_sq(t) for t in (top, second, first)).sqrt()
        assert frame_error < arb(".18")


def certify_balanced_lines(matrices, left, right):
    n = 24
    outside = [adj(left)*c*right for c in model["complement"](matrices[:6])]
    phases = [acb(0,2*arb.pi()*k/n).exp() for k in range(n)]
    for p, q in product(phases, repeat=2):
        u = acb_mat([[1], [p], [q]]) / arb(3).sqrt()
        probe = acb_mat([(adj(u)*c).entries() for c in outside])
        assert count(rows(adj(probe)*probe), arb(".39")**2) == 0
    radius = (arb(2)/3).sqrt()*arb.pi()/n
    tail = (arb(".01914")**2 + arb(".00097")**2).sqrt()
    margin = arb(".39") - radius - arb(".18") - tail
    assert margin > arb(".08")
    return margin


if __name__ == "__main__":
    for zero in (43, 195):
        matrices, _ = seven["spectral_band"](zero)
        left, right = orbit_frame(matrices)
        three_mode_bounds([adj(left)*b*right for b in matrices[:3]])
        margin = certify_balanced_lines(matrices, left, right)
        print("zero", zero, "forbidden heavy-state modulus radius >", margin.lower().str(8))
