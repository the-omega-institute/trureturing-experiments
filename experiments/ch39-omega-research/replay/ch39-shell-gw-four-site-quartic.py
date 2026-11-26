"""A four-distinct-site quartic observable in the ultraviolet singlet sector."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb

finite = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-mass-eigenvalue.py")))
q, base = finite["q"], finite["base"]


def singlet_order():
    full = base["flow"]["full"]
    s, family, _, _ = full["model"]["m"]["transported_spectra"]()
    rows, count = full["rows"], full["count"]
    positions = []
    for k, v in enumerate(np.linalg.eigvalsh(full["model"]["array"](s))):
        lo, hi = arb(float(v))-q("1e-9"), arb(float(v))+q("1e-9")
        assert count(rows(s), lo) == k and count(rows(s), hi) == k+1
        first = count(rows(family), lo-q("12e-8"))
        assert first == count(rows(family), hi+q("12e-8"))
        positions.append(k+first)
    assert positions == [0, 23, 26, 29, 42, 45]
    return positions


def certify():
    positions = singlet_order()
    events = []
    z = q("-.00675").union(q("-.0035"))
    T, B, A, margin, raw, time = base["check"](q(81, 977)**5, mass_path=events, bulk_portal=z)
    Y = (12*T).sqrt()
    heavy = q(12).sqrt()*raw.sqrt()*q(81, 977)**5*(6*T*time).exp()
    rate = 12*A+3*q(".020")+4*Y**2
    E = (8*Y*B+96*Y**3)*heavy*((rate*time).exp()-1)/rate
    assert E < min(6*margin, q(".02"))
    self_split, squarefree = [], []
    probability = q(5, 18)
    for op, *args in events:
        if op == "threshold":
            k, _, quartic, _ = args
        else:
            _, quartic, _, _, _ = args
        if k >= positions[1]:
            self_split.append((quartic["l", 0]-quartic["l", 2]-E/3).lower())
        if k >= positions[-1]:
            c = sum(quartic[t, (2, 2)] for t in ("a", "b"))
            mixed = sum(quartic[t, (0, 2)] for t in ("a", "b"))
            beta = quartic["l", 0]-mixed+c/2
            gamma = quartic["l", 2]-c/2
            assert gamma > 0
            # Each nontrivial flat singlet phase column has zero sum.
            # Its square-free coefficient is (36+|sum p_i²|²)/6^4.
            functional = probability*beta+q(5*36, 6**4)*gamma
            squarefree.append((functional-E*probability.sqrt()/6).lower())
    split = min(self_split)
    mixed = min(squarefree)
    distance = mixed*6/probability.sqrt()
    print("from heavy threshold 23: light/heavy singlet self gap above", split.str(8))
    print("from heavy threshold 45: four-site quartic above", mixed.str(8),
          "Frobenius distance from any three-site quartic above", distance.lower().str(8))
    assert split > 0 and mixed > q(".009") and distance > q(".10")


if __name__ == "__main__":
    certify()
