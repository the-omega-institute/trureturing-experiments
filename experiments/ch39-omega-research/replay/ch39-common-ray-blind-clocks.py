"""枚举三条射线首站的全部轻带盲点，并在物理钟环面比较五种荷。"""
from itertools import combinations, product
from math import gcd
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, acb_mat, arb, ctx


def first_clock_roots(fixed_row, first_row, second_row):
    # x, y are the conjugate effective phases on the target row.
    # A + B x + C y = 0 implies Re(conj(A) B x) = t.
    a, b, c = fixed_row, first_row, second_row
    if abs(a) > abs(b) + abs(c):
        return []
    d = a.conjugate() * b
    t = (abs(c)**2 - abs(a)**2 - abs(b)**2) / 2
    discriminant = abs(d)**2 - t**2
    assert abs(c) > 0 and discriminant > 0
    height = discriminant.sqrt()
    return [acb(t, sign * height) / d for sign in (-1, 1)]


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(root / "ch39-native-light-projector-pump.py"))
    ctx.prec = 160
    adj = pump["adjoint"]
    frames = list(pump["rotation_bases"](shell))
    heads = set(shell["boundary_sites"]())
    assert heads == {6, 7, 8}
    clock_edges = [(i, j) for i, j, d, _ in shell["geo"]["edges"] if d in (724, 804)]
    assert all(i not in heads and j in heads for i, j in clock_edges)
    pole = shell["soft"]["pole"]
    tuned = set(shell["soft"]["star"]["orbit"])
    assert {(min(t,pole), max(t,pole)) for t in tuned} <= {
        (i,j) for i,j,_,_ in shell["geo"]["edges"]
    }
    deficit = arb(0).union(arb("1e-8"))
    candidates = {}
    for q, active in ((1,2), (-4,1), (2,2), (-3,1), (6,2)):
        base, links = pump["clock_matrices"](shell, q)
        for t in tuned:
            base[pole,pole] -= deficit
            base[t,t] -= deficit
            base[pole,t] += deficit
            base[t,pole] += deficit

        for k, frame in enumerate(frames):
            boundary = [j for j in range(frame.ncols()) if any(
                abs(frame[i,j]) > 0 for i in heads)]
            (h,) = boundary
            interior = [j for j in range(frame.ncols()) if j != h]
            fixed = adj(frame) * base * frame
            deleted = acb_mat([[fixed[i,j] for j in interior] for i in interior])
            entries, gap = pump["first_line"](deleted)
            assert gap > arb("2.31")
            line = acb_mat([[z] for z in entries])
            assert (adj(line) * deleted * line)[0,0].real < arb("8.685")
            matrices = [fixed] + [adj(frame) * adj(link) * frame for link in links]
            coefficients = [sum(matrix[h,i] * z for i, z in zip(interior, entries, strict=True))
                            for matrix in matrices]
            roots = first_clock_roots(*coefficients)
            assert len(roots) == (2 if k == active else 0)
            if roots:
                candidates[q] = roots
                print("charge, character, physical blind points", q, k, 2*q*q)

    for q, r in combinations(candidates, 2):
        g = gcd(q, r)
        separation = min(abs(x**(r//g) - u**(q//g)).lower()
                         for x, u in product(candidates[q], candidates[r]))
        assert separation > arb("0.3001")
        print("charges, first-clock chord lower bound", q, r, separation.str(10))


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve())
