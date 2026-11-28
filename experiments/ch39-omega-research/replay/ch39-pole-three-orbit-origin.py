"""Certify the common-clock origin box after one loss from each of three orbits."""

from itertools import combinations, product
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx
import numpy as np


def numerical(matrix):
    return np.array([[complex(value) for value in row] for row in matrix.tolist()])


def dyadic(matrix):
    return acb_mat([[acb(float(value.real), float(value.imag)) for value in row]
                    for row in matrix])


def frobenius(matrix):
    return sum(abs(value)**2 for value in matrix.entries()).sqrt()


def identity(size):
    return acb_mat([[int(i == j) for j in range(size)] for i in range(size)])


def negative_pivots(matrix):
    """Count negative eigenvalues using interval Hermitian congruence."""
    lower = acb_mat(matrix)
    negative = 0
    for k in range(lower.nrows()):
        pivot = lower[k, k].real
        assert pivot < 0 or pivot > 0
        negative += int(pivot < 0)
        for i in range(k + 1, lower.nrows()):
            for j in range(k + 1, i + 1):
                lower[i, j] -= lower[i, k] * lower[j, k].conjugate() / pivot
    return negative


def trial_readout(shell, pump):
    adjoint = pump["adjoint"]
    frames = list(pump["rotation_bases"](shell))
    width = arb(1) / 10**6
    maximum_deficit = arb(1) / 10**8
    separation = arb(183) / 200
    columns = []
    squared_errors = [arb(0) for _ in frames]
    for charge in (1, -4, 2, -3, 6):
        origin, links = pump["clock_matrices"](shell, charge)
        clock_norm = arb(0)
        for link in links:
            origin += link + adjoint(link)
            row_bounds = [sum(abs(link[i, j]) + abs(link[j, i]) for j in range(19))
                          for i in range(19)]
            clock_norm += max(bound.upper() for bound in row_bounds)
        move = abs(charge) * width * clock_norm + 4 * maximum_deficit
        for sector, frame in enumerate(frames):
            frame_adjoint = adjoint(frame)
            center = frame_adjoint * origin * frame
            # Weyl transfers this count to the entire angle and deficit box.
            shift = arb("9.6") + move
            assert negative_pivots(center - shift * identity(center.nrows())) == 1
            energies, vectors = np.linalg.eigh(numerical(center))
            energy = arb(float(energies[0]))
            assert energy < arb("8.685")
            trial = dyadic(vectors[:, :1])
            lifted = frame * trial
            directions = []
            remainder_factors = []
            for link in links:
                forward = frame_adjoint * link * lifted
                backward = frame_adjoint * adjoint(link) * lifted
                directions.append(acb(0, charge) * (forward - backward))
                remainder_factors.append(frobenius(forward) + frobenius(backward))
            # The norm of a linear map of the angle pair is convex on the box.
            # The four corner norms reduce to two by simultaneous sign reversal.
            linear = width * max(frobenius(directions[0] + directions[1]).upper(),
                                 frobenius(directions[0] - directions[1]).upper())
            remainder = (charge * width)**2 * sum(remainder_factors) / 2
            size = frobenius(trial)
            residual = frobenius(center * trial - energy * trial)
            miss = (residual + linear + remainder + 4 * maximum_deficit * size) / separation
            assert miss < size
            # The radial error of the normalized orthogonal projection is
            # bounded using ||v|| - ||Pv|| = ||v-Pv||^2 / (||v|| + ||Pv||).
            radial_error = abs(size - 1) + miss**2 / size
            squared_errors[sector] += (miss + radial_error)**2
            columns.append(lifted.entries())
    # The unrestricted errors act between orthogonal cyclic sectors.
    transfer = max(bound.sqrt().upper() for bound in squared_errors)
    return acb_mat(columns).transpose(), transfer


def erasure_representatives(sites, pole):
    rotation = [sites.index((y, z, x)) for x, y, z in sites]
    assert len(sites) == 19 and rotation[pole] == pole
    orbits = [(i, rotation[i], rotation[rotation[i]]) for i in range(19)
              if i != pole and i == min(i, rotation[i], rotation[rotation[i]])]
    assert len(orbits) == 6
    cases = []
    physical_masks = set()
    for a, b, c in combinations(orbits, 3):
        for j, k in product(range(3), repeat=2):
            cases.append((a[0], b[j], c[k]))
            for turn in range(3):
                erased = frozenset((a[turn], b[(j + turn) % 3], c[(k + turn) % 3]))
                assert len(erased) == 3 and erased not in physical_masks
                physical_masks.add(erased)
    assert len(cases) == 180 and len(physical_masks) == 540
    return cases


def verify(source):
    assert __debug__
    ctx.prec = 144
    shell = run_path(str(source / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(source / "ch39-native-light-projector-pump.py"))
    sites = shell["geo"]["sites"]
    pole = shell["soft"]["pole"]
    tuned = shell["soft"]["star"]["orbit"]
    assert len(set(tuned)) == 3 and pole not in tuned
    assert all(any({i, j} == {pole, station} for i, j, _, _ in shell["geo"]["edges"])
               for station in tuned)
    readout, transfer = trial_readout(shell, pump)
    cases = erasure_representatives(sites, pole)
    minimum = arb(1)
    input_identity = identity(15)
    for erased in cases:
        keep = [i for i in range(19) if i not in erased]
        assert len(keep) == 16 and pole in keep
        retained = acb_mat([[readout[i, j] for j in range(15)] for i in keep])
        left_inverse = dyadic(np.linalg.pinv(numerical(retained)))
        defect = frobenius(input_identity - left_inverse * retained).upper()
        assert defect < 1
        center_floor = ((1 - defect) / frobenius(left_inverse).upper()).lower()
        minimum = min(minimum, center_floor)
    floor = (minimum - transfer).lower()
    assert floor > arb("1.32e-5")
    print("Certified rotation classes:", len(cases))
    print("Trial floor:", minimum.str(12))
    print("Transfer bound:", transfer.str(12))
    print("Normalized smallest singular value:", floor.str(12))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: ch39-pole-three-orbit-origin.py PINNED_REPLAY_DIRECTORY")
    verify(Path(sys.argv[1]).resolve())
