"""Uniform readout after one or two station losses in a single three-site orbit."""

from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx
import numpy as np


CELLS = 512
MAX_DEPTH = 12
THRESHOLDS = (1000, 1250, 1500)


def cover(blocks, passed, counts, full, depth=0):
    # [physical cells, characters, Bernstein y, Bernstein x, rows, charges, real/imag]
    for k, threshold in enumerate(THRESHOLDS):
        active = (passed.sum(axis=1) < 2) & ~passed[:, k]
        if active.any():
            passed[active, k] = full["inverse_test"](blocks[active, k], threshold)
    failed = passed.sum(axis=1) < 2
    counts[depth] += int((~failed).sum())
    if failed.any():
        assert depth < MAX_DEPTH
        inherited = passed[failed]
        for half in full["split_axis"](blocks[failed], 2):
            for child in full["split_axis"](half, 3):
                cover(child, inherited.copy(), counts, full, depth + 1)


def verify(root):
    full = run_path(str(Path(__file__).with_name("ch39-five-charge-full-torus.py")))
    pump = run_path(str(root / "ch39-native-light-projector-pump.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    pump["verify"](root)
    ctx.prec = 160
    scale = full["SCALE"]
    nets = []
    for k, threshold in enumerate(THRESHOLDS):
        full["verify_character"](pump, shell, k)
        net = np.empty((CELLS, CELLS, 3, 3, 6, 5, 2), dtype=np.int32)
        interpolation, off_line = full["fill_controls"](pump, shell, k, net)
        rounding = arb(50).sqrt() * (13 + 4 * MAX_DEPTH) / scale
        beta = (
            arb(threshold) / scale - interpolation - rounding - arb(5).sqrt() * off_line
        ) / arb("1.5")
        assert beta > arb(".000004")
        nets.append(net)
        print(
            "character and five-row certified floor on accepted boxes",
            k,
            beta.str(10),
            flush=True,
        )
    for omitted in range(6):
        counts = [0] * (MAX_DEPTH + 1)
        for y in range(CELLS):
            blocks = np.stack(
                [np.delete(net[y], omitted, axis=3) for net in nets], axis=1
            )
            cover(blocks, np.zeros((CELLS, 3), dtype=bool), counts, full)
        assert (
            sum(n * 4 ** (MAX_DEPTH - d) for d, n in enumerate(counts))
            == CELLS**2 * 4**MAX_DEPTH
        )
        print("omitted orbit and leaf counts by depth", omitted, counts, flush=True)
    alpha, beta = arb(".00004"), arb(".000004")
    bound = 1 / (1 / beta + (1 + arb(3).sqrt() + arb(10).sqrt() / beta) / alpha)
    assert bound > arb("5e-11")
    print("one surviving site on the damaged orbit, singular floor", bound.str(10))


if __name__ == "__main__":
    if len(sys.argv) != 2 or not __debug__:
        raise SystemExit("Run without -O and supply the pinned replay directory.")
    verify(Path(sys.argv[1]).resolve())
