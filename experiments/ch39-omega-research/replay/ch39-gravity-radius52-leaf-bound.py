"""Distance from the one-pass cutoff-52 family to the sealed cutoff-56 operators."""

import argparse
from pathlib import Path
from runpy import run_path

from flint import acb, arb


def verify(replay):
    def load(name):
        return run_path(str(replay / f"ch39-{name}.py"))

    fifth, columns, levels = load("fifth-full-positive-joint")["verify"](
        replay, return_operator="spectral")
    station = fifth*columns
    ordered = load("two-shell-ordered-bond-gates")
    sealed = load("three-shell-gravity-radius56")
    pairs = ordered["LAYERS"][0][1]
    parents, leaves = zip(*pairs, strict=True)
    assert tuple(parents) == (6, 7, 8) and tuple(leaves) == (11, 12, 18)
    assert all(site not in leaves for _, edges in ordered["LAYERS"][1:]
               for pair in edges for site in pair)
    layers = ordered["LAYERS"] + ((3, ((11, 9), (12, 15), (18, 16))),
                                   (3, ((11, 13), (12, 10), (18, 17))))
    parameters, active = sealed["PARAMETERS"], sealed["ACTIVE"]
    box = [arb(parameters[k]) + arb(0, arb("0.0002")) for k in active]
    transported = ordered["transport"](station, box, parameters=parameters,
                                       active=active, layers=layers)[0]

    # Test the top eigenvector, the second of the two heavy columns.
    test = (acb(2), acb(2, 3), acb(-1, 2))
    size = sum(abs(z)**2 for z in test).sqrt()

    def projection(matrix, sites):
        return sum(z.conjugate()*matrix[site, 1] for z, site in zip(test, sites, strict=True))

    limit = (abs(projection(station, parents))**2 + abs(projection(station, leaves))**2).sqrt()
    certified = abs(projection(transported, leaves))
    distance = (certified-limit)/size
    low, heavy, top = levels
    assert low < heavy < top
    assert distance > arb("0.10")
    # For two unit vectors the sine of their angle is d*sqrt(1-d^2/4),
    # increasing for the phase-minimized distance 0 <= d <= sqrt(2).
    d = arb("0.10")
    spectral_bound = (top-heavy)*d*(1-d**2/4).sqrt()
    assert spectral_bound > arb("0.0325")
    print("Separation of top unit eigenvectors, up to phase:", distance.lower())
    print("Induced operator-norm separation exceeds", spectral_bound.lower())


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("replay", type=Path,
                   help="replay directory of the pinned cutoff-56 experiment")
    verify(p.parse_args().replay.resolve())
