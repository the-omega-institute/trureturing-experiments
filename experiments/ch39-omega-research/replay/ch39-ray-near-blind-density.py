"""在三条射线首站的孤立盲点处计算小奇异值事件的面积系数。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb_mat, arb, ctx


def density_at_charge(shell, pump, frame, charge, clock_roots):
    adj = pump["adjoint"]
    base, links = pump["clock_matrices"](shell, charge)
    fixed = adj(frame) * base * frame
    rows = [adj(frame) * adj(link) * frame for link in links]
    n = frame.ncols()
    h, = [j for j in range(n) if abs(frame[6,j]) > 0]
    keep = [j for j in range(n) if j != h]
    deleted = acb_mat([[fixed[i,j] for j in keep] for i in keep])
    ground, gap = pump["first_line"](deleted)
    assert gap > 2
    line = acb_mat(n, 1)
    for j, value in zip(keep, ground, strict=True):
        line[j,0] = value
    energy = (adj(line) * fixed * line)[0,0].real
    a, b, c = [(op * line)[h,0] for op in (fixed, *rows)]
    identity = acb_mat([[int(i == j) for j in range(n)] for i in range(n)])
    projector = line * adj(line)
    roots = clock_roots(a, b, c)
    assert len(roots) == 2
    density = arb(0)
    for x in roots:
        y = -(a + b*x)/c
        phase = rows[0]*x + rows[1]*y
        matrix = fixed + phase + adj(phase)
        # e_h is perpendicular to the extended dark state. On its
        # orthogonal complement, the regularized inverse is the reduced resolvent.
        inverse = (matrix - identity*energy + projector).inv()
        green = inverse[h,h].real
        triangle = abs(((b*x).conjugate()*c*y).imag)
        assert green > 0 and triangle > 0
        density += 1/(4*arb.pi()*triangle*green**2)
    return density


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(root / "ch39-native-light-projector-pump.py"))
    prior = run_path(str(Path(__file__).with_name("ch39-common-ray-blind-clocks.py")))
    ctx.prec = 160
    frames = list(pump["rotation_bases"](shell))
    cases = ((1,2,"25.99278856"), (-4,1,"17.37426488"), (2,2,"26.16959086"),
             (-3,1,"12.59217795"), (6,2,"14.28551665"))
    total = arb(0)
    for charge, character, reference in cases:
        density = density_at_charge(shell, pump, frames[character], charge, prior["first_clock_roots"])
        assert abs(density - arb(reference)) < arb("1e-8")
        print("荷、近盲频率二次项系数", charge, density.str(13))
        total += density
    assert abs(total - arb("96.4143389021")) < arb("1e-10")
    print("至少一种荷近盲的二次项系数", total.str(14))


if __name__ == "__main__":
    verify(Path(sys.argv[1]).resolve())
