"""Bound the least singular value of every fifteen-site readout at one zeta zero."""

from itertools import combinations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


CHARGES = (1, -4, 2, -3, 6)
ZERO_INDEX = 2_306_888_100_085


def readout_columns(shell, pump, height):
    """Use the physical Laplacian for the softened model and the clock pencil for the flat one."""
    frames = list(pump["rotation_bases"](shell))
    phases = shell["geo"]["phases"](height)
    flat_columns, actual_columns = [], []
    adjoint = pump["adjoint"]
    for charge in CHARGES:
        base, links = pump["clock_matrices"](shell, charge)
        flat = base
        for label, link in zip((724, 804), links, strict=True):
            term = link * phases[label]**charge
            flat += term + adjoint(term)
        actual = shell["laplacian"](charge, height)
        for frame in frames:
            for operator, columns in ((flat, flat_columns), (actual, actual_columns)):
                line, gap = pump["first_line"](adjoint(frame) * operator * frame)
                assert gap > 1
                unit = frame * acb_mat([[z] for z in line])
                columns.append(unit.entries())
    return flat_columns, actual_columns


def frobenius_upper(matrix):
    return sum(z.abs_upper()**2 for z in matrix.entries()).sqrt().upper()


def certify(columns):
    identity = acb_mat([[int(i == j) for j in range(15)] for i in range(15)])
    smallest = arb(1)
    checked = 0
    for sites in combinations(range(19), 15):
        matrix = acb_mat([[column[i] for column in columns] for i in sites])
        inverse = matrix.mid().inv().mid()
        floor = ((1 - frobenius_upper(identity - inverse * matrix))
                 / frobenius_upper(inverse)).lower()
        assert floor > arb("0.00000311736"), sites
        smallest = min(smallest, floor)
        checked += 1
    assert checked == 3876
    return checked, smallest


def main(inputs):
    shell = run_path(str(inputs / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(inputs / "ch39-native-light-projector-pump.py"))
    ctx.prec = 192
    height = acb.zeta_zero(ZERO_INDEX).imag
    for model, columns in zip(("flat", "actual"), readout_columns(shell, pump, height), strict=True):
        count, floor = certify(columns)
        print(model, count, floor.str(14), flush=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
