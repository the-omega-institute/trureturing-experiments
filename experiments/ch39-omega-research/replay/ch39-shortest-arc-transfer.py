"""Finite transfer and cutoff exclusion for the sphere of squared radius 867."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from sympy import Matrix, eye, zeros


def verify():
    points = [p for p in product(range(30), repeat=3) if sum(x*x for x in p) == 867]
    edges = [(sum((a-b)**2 for a, b in zip(p, q)), i, j)
             for (i, p), (j, q) in combinations(enumerate(points), 2)]
    cuts = sorted({d for d, i, j in edges})
    assert cuts[:4] == [12, 32, 52, 56]
    cycle_data = []
    for cut in cuts[:4]:
        reps = list(range(19))
        degrees = Counter()
        for d, i, j in edges:
            if d <= cut:
                a, b = reps[i], reps[j]
                reps = [b if r == a else r for r in reps]
                degrees.update([i, j])
        count = sum(degrees.values())//2
        cycle_data.append((count, count-19+len(set(reps))))
        if cut >= 32:
            assert sum(v*v for v in degrees.values()) >= 12
    assert cycle_data == [(3, 0), (6, 0), (12, 0), (18, 6)]
    a = zeros(19)
    for d, i, j in edges:
        if d == 12:
            a[i, j] = a[j, i] = 1
    assert a**3 == a and (a**2).trace() == 6
    t = F(1, 19)
    epsilon = 6*t**4/(1-t**6)
    # In the forward sheet matrix, the (1,0) block is t A (I-t^3 A^3)^(-1).
    forward = (t*a*(eye(19)-t**3*a**3).inv()).trace()
    assert forward == epsilon
    # In the reciprocal case the two active vertex eigenvalues are +1 and -1.
    shift = Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
    reciprocal = sum(3*(eye(3)-s*t*(shift+shift.T)).inv()[1, 0] for s in (1,-1))
    assert reciprocal == 6*t*t/((1-t*t)*(1-4*t*t))
    upper = F(49, 10**6)
    assert F(43, 10**6) < epsilon < upper < 12*t**4/(1+(18*t)**3)
    assert reciprocal > F(1, 100)
    study = run_path(str(Path(__file__).with_name("ch39-shell-zero-mass-hierarchy.py")))
    scalar, _, fit = study["parameters"]()
    low, v, log_x, g5, g10 = fit(scalar(epsilon), scalar(1))
    assert scalar(".631") < low < scalar(".633") and scalar("7.469e14") < v < scalar("7.471e14")
    print("Ranks of the cycle spaces at squared cutoffs 12,32,52,56:", cycle_data)
    print("Forward projected response", epsilon, "reciprocal response", reciprocal)


if __name__ == "__main__":
    from flint import ctx
    ctx.prec = 160
    verify()
