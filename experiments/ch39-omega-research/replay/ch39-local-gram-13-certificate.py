"""Certify an exact rank-thirteen Gram completion at squared chord cutoff 56.

The seed table ch39-local-gram-13.txt holds a 19 x 13 integer matrix; dividing by
10**17 gives a rational matrix Y0.  A Newton-Kantorovich box around Y0 contains an
exact solution of the 153 far-pair equations, so the cutoff-56 completion W has rank
at most 13 and the three clocks are realised with 19 + 13 = 32 auxiliary states.
"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx, fmpq

root = Path(__file__).parent
prev = run_path(str(root/"ch39-local-memory-rank.py"))


def largest(values):
    """Maximum of exact rationals, or a ball containing the maximum of several balls."""
    values = list(values)
    best = values[0]
    for value in values[1:]:
        best = best.max(value) if isinstance(best, arb) else max(best, value)
    return best


def row_norm(matrix):
    return largest(sum(abs(a) for a in row) for row in matrix.tolist())


def verify():
    ctx.prec = 192
    _, _, _, g0, g1, g2 = prev["MEM"]["coefficients"]()
    delta = g2-g1*g0.solve(g1)
    points = sorted(prev["LOCAL"]["shell_grid"](1)|{(17,17,17)})
    pairs = [(i,j) for i,j in combinations(range(19), 2)
             if prev["LOCAL"]["squared"](points[i],points[j]) > 56]
    larger_scale = {(11,11,25), (13,13,23), (17,17,17)}
    scales = [8192 if tuple(sorted(p)) in larger_scale else 2048 for p in points]
    seed = arb_mat([[arb(int(a))/10**17 for a in line.split()] for line in
                    (root/"ch39-local-gram-13.txt").read_text().splitlines()])
    n, rank = seed.nrows(), seed.ncols()
    assert (n, rank, len(pairs)) == (19, 13, 153)
    near = sorted({(i,j) for i,j in combinations(range(19), 2)} - set(pairs))
    assert len(near) == 18 and 19+len(near) == 37
    assert near == [(0,2),(1,3),(4,6),(4,8),(5,7),(5,9),(6,8),(6,10),(7,9),(7,12),(8,10),
                    (9,12),(13,14),(13,16),(14,15),(14,16),(15,16),(17,18)]
    # Parameter count of a 19 x r Gram factor modulo O(r): first reaches 153 at r = 11.
    assert 19*10-45 < 153 <= 19*11-55
    # Jacobian of f(Y)_{ij} = Y_i.Y_j - s_i s_j Delta_ij over the 153 far pairs.
    jac = arb_mat(len(pairs), n*rank)
    gram = seed*seed.transpose()
    error = arb(0)
    for t,(i,j) in enumerate(pairs):
        for k in range(rank):
            jac[t,i*rank+k] = seed[j,k]
            jac[t,j*rank+k] = seed[i,k]
        error = error.max(abs(gram[i,j]-arb(delta[i,j]*scales[i]*scales[j])))
    assert error < arb("2.3e-17")
    right = jac.transpose()*(jac*jac.transpose()).inv()
    bound = row_norm(right)
    assert arb("7085.41298") < bound < arb("7085.41299") and bound < 7200
    # Fixed point u = -f(Y0) - Q(Pu) on the sup-norm ball of radius eta.
    radius = arb("1e-12")
    assert error + rank*(bound*radius)**2 < radius
    assert 2*rank*bound**2*radius < 1
    enclosed = arb_mat([[(seed[i,k]+arb(0,bound*radius))/scales[i] for k in range(rank)]
                        for i in range(n)])
    w = enclosed*enclosed.transpose()
    norm = row_norm(w)
    inverse_norm = row_norm(g0.inv())
    assert arb("1.402894e-6") < norm < arb("1.402896e-6") and inverse_norm < 315000
    assert norm*315000 < (1-arb(fmpq(98,1000)))**2
    print("Gram rank at most", rank, "right inverse norm", bound.str(8), "residual", error.str(8))
    assert arb("0.441911") < 315000*norm < arb("0.441913") and (1-arb(fmpq(98,1000)))**2 > arb("0.8136")
    print("Row norm bound", norm.str(10), "contraction product", (315000*norm).str(6),
          "auxiliary dimension at most", n+rank)
    print("Cutoff 56 interval for the auxiliary dimension:", 29, 32)


if __name__ == "__main__":
    verify()
