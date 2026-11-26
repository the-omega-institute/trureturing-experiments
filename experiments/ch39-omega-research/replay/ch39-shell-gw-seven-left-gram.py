"""Identify the first scalar band realizing every prescribed left mass Gram matrix."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

model = run_path(str(Path(__file__).with_name("ch39-shell-gw-195-flavor.py")))
bands = model["bands"]
rows = bands["m"]["rows"]
count = model["count"]


def spectral_band(index):
    if index == 195:
        single, family, vertex = model["scalar_models"]()
        bounds = ("1.1071473", "1.1071475")
        model["cuts"][7] = bounds[1]
        matrices = model["image_basis"](vertex, model["sources"](single, family, 7))
    else:
        single, family, _, _ = bands["five"]["m"]["transported_spectra"]()
        bounds = ("1.0998366", "1.0998368")
        bands["cuts"][7] = bounds[1]
        matrices = bands["orb"]["orthonormalize"](bands["images"](7))
    for cut, expected in zip(map(arb, bounds), (5, 6)):
        assert count(rows(single), cut) == 1
        assert count(rows(family), cut) == expected
    return matrices, bounds


def check(index):
    matrices, bounds = spectral_band(index)
    outside = model["complement"](matrices[:6])
    left = acb_mat([[1], [acb(0, -1)], [-1]]) / arb(3).sqrt()
    # Rows are u* C_j. Their smallest singular value bounds the norm of
    # projection of every rank-one target u v* onto the missing three modes.
    witness = acb_mat([(model["adj"](left) * c).entries() for c in outside])
    lower = arb(".51")
    assert count(rows(model["adj"](witness) * witness), lower**2) == 0
    print(f"zero {index}, seventh threshold inside {bounds}, six-mode rank-one distance > 0.51")


if __name__ == "__main__":
    for index in (43, 195):
        check(index)
