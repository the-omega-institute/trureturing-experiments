"""Bound the rational degree available to any fixed linear feedback on the charged shell."""
from collections import Counter
from runpy import run_path
from pathlib import Path
from flint import acb, arb, ctx, fmpq, fmpq_mat


def verify():
    ctx.prec = 128
    shell = run_path(str(Path(__file__).with_name("ch39-rotor-shell-returns.py")))
    nodes, adj, _ = shell["charged_shell"]()
    charges = Counter(q for row in adj for _,_,q in row)
    assert charges == {-1:6, 0:564, 1:6}
    plus = fmpq_mat(len(nodes), len(nodes))
    for i,row in enumerate(adj):
        for j,a,q in row:
            if q == 1:
                plus[i,j] = fmpq(a.numerator, a.denominator)
    native_degree = plus.rank()+plus.transpose().rank()
    assert native_degree == 12

    first, last = 10905, 12054
    zeros = acb.zeta_zeros(first,last-first+1)
    log_two, root_two = arb(2).log(), arb(2).sqrt()
    scores, signs = [], []
    for zero in zeros:
        u = root_two*(zero.imag*log_two).cos()
        if abs(u) > arb(".01"):
            signs.append(u > 0)
        scores.append((3-2*u)*(2*(113+72*u)).sqrt()/(13*(zero.imag**2+arb(1)/4)))
    flips = sum(a != b for a,b in zip(signs,signs[1:]))
    assert flips == 213
    assert all(not a.overlaps(b) for a,b in zip(scores,scores[1:]))
    slopes = [b>a for a,b in zip(scores,scores[1:])]
    extrema = sum(a!=b for a,b in zip(slopes,slopes[1:]))
    assert extrema >= 425 and sum(a and not b for a,b in zip(slopes,slopes[1:])) == 213
    assert 2*106 < flips and 4*106-2 < extrema
    print("Intrinsic charged ranks 6+6; dyadic sign flips", flips,
          "and adjacent-score direction reversals", extrema, "on", len(zeros), "guaranteed labels")


if __name__ == "__main__":
    verify()
