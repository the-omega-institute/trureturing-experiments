"""Odd powers on an opposite pair at Q force longer bilinear currents for Gaussian exchange."""
from pathlib import Path
from runpy import run_path

from flint import arb

finite = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-mass-eigenvalue.py")))
q, base = finite["q"], finite["base"]


def certify():
    events = []
    z = q("-.00675").union(q("-.0035"))
    T, B, A, margin, raw, time = base["check"](q(81, 977)**5, mass_path=events, bulk_portal=z)
    Y = (12*T).sqrt()
    heavy = q(12).sqrt()*raw.sqrt()*q(81, 977)**5*(6*T*time).exp()
    rate = 12*A+3*q(".020")+4*Y**2
    E = (8*Y*B+96*Y**3)*heavy*((rate*time).exp()-1)/rate
    assert E < min(6*margin, q(".02"))
    signal = q(1)
    for op, *args in events:
        if op == "threshold":
            k, groups, quartic, _ = args
        else:
            groups, quartic, _, _, _ = args
        if k == 53:
            assert sum(groups.values()) == 54
            old = quartic["a", (0, 1)]-quartic["a", (1, 2)]
            bulk = quartic["a", (0, 2)]-quartic["a", (2, 2)]
            assert old > bulk
            signal = arb.min(signal, ((old+5*bulk)/18-E/q(3).sqrt()).lower())
    # Real coordinates, J=svec(phi phi^T), A4=Z/24.
    gram_lower = q(3, 6)-(2*A+E)/24
    assert gram_lower > q(".12")
    print("opposite-pair t*r*h² coefficient exceeds", signal.lower().str(8),
          "Gaussian SOS completion Gram floor for C=3 exceeds", gram_lower.lower().str(8))
    assert signal > q(".024")


if __name__ == "__main__":
    certify()
