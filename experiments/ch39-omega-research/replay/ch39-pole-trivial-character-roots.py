"""Verify local roots of the trivial-character minor and full pole readout rank."""
from itertools import permutations
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, arb_mat, ctx
import numpy as np

CHARGES = (1, -4, 2, -3, 6)
CASES = (
    ((0, 2), '6.17447975051', '0.466046496704'),
    ((0, 5), '3.23224491967', '2.94247163631'),
    ((1, 2), '0.108705556666', '5.81713881048'),
    ((1, 5), '3.05094038751', '3.34071367087'),
    ((2, 3), '0.20834108299', '6.18755047017'),
    ((2, 4), '6.07484422419', '0.0956348370158'),
    ((2, 5), '0.187503289796', '5.61079842748'),
    ((3, 4), '0.00764035285244', '3.14798884339'),
    ((3, 5), '3.15132980725', '6.18823521815'),
    ((4, 5), '3.13185549993', '0.0949500890266'),
)


def infinity_norm(matrix):
    return max(sum(entry.abs_upper() for entry in row) for row in matrix.tolist())


def frobenius(matrix):
    return sum(abs(entry)**2 for entry in matrix.entries()).sqrt()


def numeric(matrix):
    return np.array([[complex(entry) for entry in row] for row in matrix.tolist()])


def dyadic(matrix):
    return acb_mat([[acb(float(z.real), float(z.imag)) for z in row] for row in matrix])


def verify(source, verifier):
    ctx.prec = 200
    common = run_path(str(verifier / 'ch39-common-charge-amplitudes.py'))
    shell = run_path(str(source / 'ch39-three-ray-chiral-boundary.py'))
    pump = run_path(str(source / 'ch39-native-light-projector-pump.py'))
    adjoint, first_line = pump['adjoint'], pump['first_line']
    frames = list(pump['rotation_bases'](shell))
    sites, pole = shell['geo']['sites'], shell['soft']['pole']
    rotation = [sites.index((b, c, a)) for a, b, c in sites]
    orbits = [(i, rotation[i], rotation[rotation[i]]) for i in range(19)
              if i != pole and i == min(i, rotation[i], rotation[rotation[i]])]
    assert [sites[group[0]] for group in orbits] == [
        (1, 5, 29), (1, 29, 5), (13, 13, 23), (7, 17, 23), (7, 23, 17), (11, 11, 25)]
    trivial = frames[0]
    orbit_rows = [next(j for j in range(7) if not trivial[group[0], j].is_zero()) for group in orbits]
    pole_row = next(j for j in range(7) if not trivial[pole, j].is_zero())
    assert len(set(orbit_rows + [pole_row])) == 7
    star = acb_mat(19, 19)
    for leaf in shell['soft']['star']['orbit']:
        star[pole, pole] -= 1
        star[leaf, leaf] -= 1
        star[pole, leaf] += 1
        star[leaf, pole] += 1
    drift = adjoint(trivial) * star * trivial
    operators = [pump['clock_matrices'](shell, charge) for charge in CHARGES]
    for _, links in operators:
        for link in links:
            symmetric = link + adjoint(link)
            assert max(sum(not symmetric[i, j].is_zero() for j in range(19)) for i in range(19)) == 2
    radius, maximum_deficit = arb('3e-8'), arb('1e-8')
    deficit_interval = arb(0).union(maximum_deficit)
    gap = arb(183) / 200
    eye2 = arb_mat([[1, 0], [0, 1]])
    eye15 = acb_mat([[int(i == j) for j in range(15)] for i in range(15)])

    def system(angles, deficit, kept, pivots=None):
        columns, chosen = [], []
        for i, (charge, (base, links)) in enumerate(zip(CHARGES, operators, strict=True)):
            iq = acb(0, charge)
            a, b = [adjoint(trivial) * link * trivial * (iq*angle).exp()
                    for angle, link in zip(angles, links, strict=True)]
            h = adjoint(trivial) * base * trivial + a + adjoint(a) + b + adjoint(b) + deficit*drift
            if pivots is None:
                line = first_line(h)[0]
                pivot = max(range(7), key=lambda j: float(abs(line[j])))
            else:
                pivot = pivots[i]
            slopes = [(a-adjoint(a))*iq, (b-adjoint(b))*iq, drift]
            columns.append(common['eigenline_derivatives'](h, slopes, first_line, pivot))
            chosen.append(pivot)
        f, (dx, dy, ds) = common['determinant_derivatives'](columns, kept)
        return (arb_mat([[f.real], [f.imag]]),
                arb_mat([[dx.real, dy.real], [dx.imag, dy.imag]]),
                arb_mat([[ds.real], [ds.imag]]), chosen)

    for pair, x, y in CASES:
        center = [arb(x), arb(y)]
        box = [angle + arb(0, radius.upper()) for angle in center]
        kept = [row for i, row in enumerate(orbit_rows) if i not in pair] + [pole_row]
        f, j, _, pivots = system(center, arb(0), kept)
        inverse = j.mid().inv()
        assert not inverse.det().contains(0)
        _, derivative, drift_bound, _ = system(box, deficit_interval, kept, pivots)
        ratio = (infinity_norm(inverse*f) + maximum_deficit*infinity_norm(inverse*drift_bound))/radius
        ratio += infinity_norm(eye2 - inverse*derivative)
        assert ratio < arb('.4'), (pair, ratio)

        columns, errors = [], []
        for charge, (base, links) in zip(CHARGES, operators, strict=True):
            h = base
            for angle, link in zip(center, links, strict=True):
                part = acb(0, charge*angle).exp() * link
                h = h + part + adjoint(part)
            for frame, soft_norm in zip(frames, (4, 1, 1), strict=True):
                block = adjoint(frame) * h * frame
                values, vectors = np.linalg.eigh(numeric(block))
                v = dyadic(vectors[:, :1])
                energy, size = arb(float(values[0])), frobenius(v)
                assert energy < arb('8.685')
                move = 4*abs(charge)*radius + soft_norm*maximum_deficit
                off_line = (frobenius(block*v - v*energy) + move*size)/gap
                assert size > off_line
                errors.append(2*off_line + abs(size - 1))
                columns.append((frame*v).entries())
        transfer = sum(error**2 for error in errors).sqrt().upper()
        floor = arb(1)
        for doubled, singleton in permutations(pair):
            for site in orbits[singleton]:
                erased = {*orbits[doubled][1:], site}
                readout = acb_mat([[column[i] for column in columns] for i in range(19) if i not in erased])
                left = dyadic(np.linalg.pinv(numeric(readout)))
                defect = frobenius(eye15 - left*readout).upper()
                assert defect < 1
                floor = min(floor, ((1-defect)/frobenius(left).upper()).lower() - transfer)
        assert floor > arb('6.61e-5'), (pair, floor)
        print(pair, 'self-map ratio', ratio.str(8), 'normalized lower bound', floor.str(10))


if __name__ == '__main__':
    assert __debug__
    verify(*(Path(arg).resolve() for arg in sys.argv[1:]))
