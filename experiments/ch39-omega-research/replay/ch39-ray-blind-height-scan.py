import argparse
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

parser = argparse.ArgumentParser(description="Certify finite-height exclusion at the ray readout.")
parser.add_argument("replay_dir", type=Path, help="replay directory at d0c5ab50ebaf763fea12e2629cb6760205448f08")
parser.add_argument("--grid-bits", type=int, default=180)
args = parser.parse_args()
if not __debug__:
    parser.error("run without -O because the verification uses assertions")
if not 128 <= args.grid_bits <= 220:
    parser.error("--grid-bits must be between 128 and 220")
root = args.replay_dir.resolve()
shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
pump = run_path(str(root / "ch39-native-light-projector-pump.py"))
ctx.prec = 256

adjoint = pump["adjoint"]
frames = list(pump["rotation_bases"](shell))
heads = set(shell["boundary_sites"]())
assert heads == {6, 7, 8}
clock_edges = [(i, j) for i, j, distance, _ in shell["geo"]["edges"] if distance in (724, 804)]
assert all(i not in heads and j in heads for i, j in clock_edges)
pole = shell["soft"]["pole"]
tuned = set(shell["soft"]["star"]["orbit"])
tau = 2 * arb.pi()
theta = shell["geo"]["theta"]
r1 = (theta(724) - theta(244)) / tau
r2 = (theta(804) - theta(244)) / tau
ratio = r2 / r1
charge_sector = ((1, 2), (-4, 1), (2, 2), (-3, 1), (6, 2))
modulus = 2**args.grid_bits
ratio_numerator = int((modulus * ratio.mid()).floor().unique_fmpz())
ratio_rational = arb(ratio_numerator) / modulus


def ceil_bracket(value):
    return int(value.upper().ceil().unique_fmpz())


def floor_sum(count, modulus, slope, offset):
    total = 0
    while count:
        slope_quotient, slope = divmod(slope, modulus)
        offset_quotient, offset = divmod(offset, modulus)
        total += slope_quotient * count * (count - 1) // 2 + offset_quotient * count
        height = slope * count + offset
        if height < modulus:
            break
        count, offset = divmod(height, modulus)
        slope, modulus = modulus, slope
    return total


def integer_residue_hits(start, end, offset, width):
    count = end - start + 1
    if count <= 0:
        return 0
    shifted = offset + ratio_numerator * start

    def below(limit):
        plain = floor_sum(count, modulus, ratio_numerator, shifted)
        shifted_sum = floor_sum(count, modulus, ratio_numerator, shifted + modulus - limit)
        return count + plain - shifted_sum

    if width * 2 + 1 >= modulus:
        return count
    return below(width + 1) + count - below(modulus - width)


def enclose_roots(charge, sector, U):
    base, links = pump["clock_matrices"](shell, charge)
    height = arb(U).union(arb(2 * U))
    offset = shell["soft"]["point"]["study"]["offset"]
    deficit = arb(23) / 25 * arb(2)**-20 / (height**2 + offset)
    assert 0 < deficit < arb("1e-8")
    for t in tuned:
        base[pole, pole] -= deficit
        base[t, t] -= deficit
        base[pole, t] += deficit
        base[t, pole] += deficit
    frame = frames[sector]
    (boundary,) = [j for j in range(frame.ncols()) if any(abs(frame[i,j]) > 0 for i in heads)]
    interior = [j for j in range(frame.ncols()) if j != boundary]
    fixed = adjoint(frame) * base * frame
    deleted = acb_mat([[fixed[i, j] for j in interior] for i in interior])
    entries, gap = pump["first_line"](deleted)
    assert gap > arb("2.31")
    line = acb_mat([[entry] for entry in entries])
    assert (adjoint(line) * deleted * line)[0,0].real < arb('8.685')
    matrices = [fixed] + [adjoint(frame) * adjoint(link) * frame for link in links]
    a, b, c = [sum(matrix[boundary, i] * z for i, z in zip(interior, entries, strict=True))
               for matrix in matrices]
    if abs(a) > abs(b) + abs(c):
        return
    cross = a.conjugate() * b
    projection = (abs(c)**2 - abs(a)**2 - abs(b)**2) / 2
    disc = abs(cross)**2 - projection**2
    assert abs(c) > 0 and disc > 0
    for sign in (-1, 1):
        x = acb(projection, sign * disc.sqrt()) / cross
        y = -(a + b*x)/c
        assert abs(x).contains(1) and abs(y).contains(1)
        if charge > 0:
            x, y = x.conjugate(), y.conjugate()
        yield x.arg() / tau, y.arg() / tau


def classify(U, q, k, active):
    height = arb(U).union(arb(2 * U))
    roots = list(enclose_roots(q, k, U))
    assert len(roots) == (2 if active else 0)
    for branch, (x, y) in enumerate(roots):
        range_ball = abs(q) * r1 * height - x
        start = int(range_ball.lower().ceil().unique_fmpz())
        end = int(range_ball.upper().floor().unique_fmpz())
        if start > end:
            continue
        constant = ratio * x - y
        numerator = int((modulus * constant.mid()).floor().unique_fmpz())
        error = abs(constant - arb(numerator)/modulus) + abs(ratio - ratio_rational)*max(abs(start), abs(end))
        width = ceil_bracket(modulus * error)
        hits = integer_residue_hits(start, end, numerator, width)
        assert hits == 0, (U, q, k, branch, start, end, hits)


assert acb.zeta_zero(1).imag > 14
for exponent in range(80):
    lower_height = 14 * 2**exponent
    for charge, active_sector in charge_sector:
        for sector in range(3):
            classify(lower_height, charge, sector, sector == active_sector)
    print("scanned", lower_height, 2 * lower_height, flush=True)
print("No continuous-height hits", 14, 14 * 2**80, "grid bits", args.grid_bits)
