"""Two-sided certified strength for a free common field on the physical quark orbit."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb

common = run_path(str(Path(__file__).with_name("ch39-shell-gw-195-common-field.py")))
adj, hs, rows, count = (common[key] for key in ("adj", "hs", "rows", "count"))
orbit = common["flavor"]["bands"]["orb"]
seeds = {
    43: (
        """-78983542 43518684 -33850889 26855601 687395 -381746
-332699 -262064 -1008231 782300 -99931206 3456326
-42464407 8019914 84586418 -31247453 -993272 -450618""",
        """-42221518 -73883016 -846551 -402138 41909300 -31644519
-42216410 -31239756 291643 -497528 -41513112 74284462
-663400 235778 84572666 -53350288 369210 -760885""",
        """-51835080 -6817728 -101837 869699 -67741141 -51739870
2949602 -85188236 166300 -1043208 -27606387 44396349
-123737 -944508 48200393 87606161 803828 573708""",
    ),
    195: (
        """-77642765 42341609 -36709172 28769285 -385868 1845323
43649499 16365319 -86804515 -16885597 -371926 -2571382
2129066 1074172 -352593 -2119755 -62664772 77864125""",
        """-42343297 -73974147 31343103 41740746 2761878 1600900
-41978640 -31093286 -74227936 -41938478 -914164 1242257
1134868 2141092 354508 -2563567 65936004 75099092""",
        """-27799111 2519536 94603650 16150992 2633534 -1794650
8248410 -95638350 9498424 -26307275 -131699 1720336
106472 2473799 1557147 -2139496 22322071 97409404""",
    ),
}


def frames(index):
    result = []
    for seed in seeds[index]:
        integers = list(map(int, seed.split()))
        # A common positive denominator drops out of Gram-Schmidt.
        columns = []
        for j in range(3):
            v = acb_mat([[acb(*integers[6*i+2*j:6*i+2*j+2])] for i in range(3)])
            for col in columns:
                v -= hs(col, v)*col
            norm2 = hs(v, v).real
            assert norm2 > 0
            columns.append(v/norm2.sqrt())
        result.append(acb_mat([[v[i, 0] for v in columns] for i in range(3)]))
    return result


def readouts(index):
    F = common["old"]["F"] if index == 43 else common["operators"]()[1]
    phase = common["native"]["relative_phase"](arb(676))
    return F, common["conjugate"](F, phase)


def limits(index, F, D, up, down):
    bound = arb(".010790772074" if index == 43 else ".010772466475")
    assert count(rows(adj(F)*F-60*adj(D)*D), bound) == 54
    lower = ((hs(up, up).real-60*hs(down, down).real)/bound).sqrt()
    U, A, B = frames(index)
    target = acb_mat([[z] for z in (U*up*A).entries() + (U*down*B).conjugate().entries()])
    P = acb_mat(rows(F)+rows(D))
    field = adj(P)*(P*adj(P)).solve(target)
    strength = hs(field, field).real.sqrt()
    return lower, strength, field


if __name__ == "__main__":
    V, u, d = orbit["reference"]()
    for n, lo, hi in ((43,"8.2282703","8.2284765"),(195,"8.2352584","8.2356297")):
        lower, upper, field = limits(n, *readouts(n), u*arb(".8616"), V*d*arb(".01401"))
        assert arb(lo) < lower and upper < arb(hi)
        print("zero", n, "infimum bounded below by", lower.str(13), "witness norm", upper.str(13))
