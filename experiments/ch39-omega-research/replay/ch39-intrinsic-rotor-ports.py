"""Intrinsic maximal charged-return ports, their polynomial selector and weighting threshold."""
from collections import Counter
from fractions import Fraction
from math import prod
from pathlib import Path
from runpy import run_path
from flint import acb, arb, ctx


def ball(x):
    return arb(x.numerator)/x.denominator


def verify():
    ctx.prec = 192
    shell = run_path(str(Path(__file__).with_name("ch39-rotor-shell-returns.py")))
    nodes, adj, root = shell["charged_shell"]()
    rows = []
    for v in range(len(nodes)):
        _, _, three, four = shell["returns"](adj, v)
        rows.append((three.get(1, Fraction()), four.get(1, Fraction()), four.get(2, Fraction())))
    levels = sorted(set(b for b, _, _ in rows))
    high = levels[-1]
    assert len(levels) == 5 and high == rows[root][0]
    selected = {v for v, (b, _, _) in zip(nodes, rows) if b == high}
    assert selected == {(33,33,75), (33,75,33), (75,33,33)}
    for b, _, _ in rows:
        assert prod((b-s)/(high-s) for s in levels[:-1]) == int(b == high)

    parameter = Fraction(7460124,259581545)
    x = [1-1/(2*(z.imag**2+arb(1)/4)) for z in acb.zeta_zeros(10874,3)]
    grouped = Counter(rows)
    groups = []
    for (b,a,c), count in grouped.items():
        if b:
            slopes = [count*(ball(a-parameter*b)+2*ball(c)*(u+v)) for u,v in zip(x,x[1:])]
            groups.append((ball(b/high), slopes))
    groups.sort(key=lambda g:float(g[0]),reverse=True)
    assert len(groups) == 4 and groups[0][0] == 1
    assert groups[0][1][0] < 0 < groups[0][1][1]
    assert all(s < 0 for s in groups[1][1])
    assert all(s > 0 for _,pair in groups[2:] for s in pair)
    for power in range(1,215):
        right = sum(r**power*slopes[1] for r,slopes in groups)
        assert right < 0
    # At every larger power the two positive error terms shrink, as does the negative one.
    left_upper = groups[0][1][0] + sum(r**215*slopes[0] for r,slopes in groups[2:])
    right_lower = groups[0][1][1] + groups[1][0]**215*groups[1][1][1]
    assert left_upper < 0 < right_lower
    # Power zero uses every port, including those with b=0.
    assert sum(ball(a-parameter*b)+2*ball(c)*(x[0]+x[1]) for b,a,c in rows) > 0
    print("5 return amplitudes, degree-4 exact selector, 3 ports, sharp positive integer power 215")


if __name__ == "__main__":
    verify()
