"""Interval bounds on blind-line site amplitudes and the head-deleted resolvent."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path
import sys

from flint import acb_mat, arb, ctx


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(root / "ch39-native-light-projector-pump.py"))
    adjoint = pump["adjoint"]
    ctx.prec = 160
    frames = list(pump["rotation_bases"](shell))
    heads = set(shell["boundary_sites"]())
    pole = shell["soft"]["pole"]
    assert heads == {6, 7, 8} and pole == 14
    extras = sorted(set(range(19)) - heads - {pole})
    links_read_at_heads = [(i, j) for i, j, d, _ in shell["geo"]["edges"] if d in (724, 804)]
    assert all(i not in heads and j in heads for i, j in links_read_at_heads)
    deficit = arb(0).union(arb("1e-8"))
    site_floor = {i: arb(1) for i in extras}
    for charge, critical in ((1, 2), (-4, 1), (2, 2), (-3, 1), (6, 2)):
        base, links = pump["clock_matrices"](shell, charge)
        for site in shell["soft"]["star"]["orbit"]:
            base[pole, pole] -= deficit
            base[site, site] -= deficit
            base[pole, site] += deficit
            base[site, pole] += deficit
        for k, frame in enumerate(frames):
            h, = [j for j in range(frame.ncols()) if any(abs(frame[i, j]) > 0 for i in heads)]
            interior = [j for j in range(frame.ncols()) if j != h]
            fixed = adjoint(frame) * base * frame
            deleted = acb_mat([[fixed[i, j] for j in interior] for i in interior])
            entries, gap = pump["first_line"](deleted)
            assert gap > arb("2.31") and 0 < fixed[h, h].real < 20
            row_blocks = [fixed] + [adjoint(frame) * adjoint(m) * frame for m in links]
            row_norm = sum(sum(abs(m[h, i])**2 for i in interior).sqrt() for m in row_blocks)
            assert row_norm < 9
            if k == critical:
                # A site fixed by the rotation has zero coordinates in both nontrivial characters.
                assert all(frame[pole, j].is_zero() for j in range(frame.ncols()))
                for i in extras:
                    amplitude = sum(frame[i, j] * z for j, z in zip(interior, entries, strict=True))
                    site_floor[i] = min(site_floor[i], abs(amplitude).lower())
            else:
                a, b, c = [sum(m[h, j] * z for j, z in zip(interior, entries, strict=True))
                           for m in row_blocks]
                assert abs(a) - abs(b) - abs(c) > arb("0.149")
    assert all(z > arb("0.0174") for z in site_floor.values())
    print("Blind-state lower bounds per extra site:",
          {i: int((v*10**6).floor().unique_fmpz()) for i, v in site_floor.items()})

    # Exact constants in the proof, with head Fourier coordinates normalized to unit norm.
    b, small = Fraction(1, 400), Fraction(1, 20000)
    overlap, amplitude = Fraction(999, 1000), Fraction(17, 1000)
    floor = Fraction(4, 100000)
    assert 1 - 17*b*b > overlap**2
    assert Fraction(149, 1000)*overlap > 56*b
    assert Fraction(174, 10000)*overlap - 4*small > amplitude
    assert floor < small and floor**2 * (amplitude**2 + (1+b)**2) < (amplitude*b)**2


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve())
