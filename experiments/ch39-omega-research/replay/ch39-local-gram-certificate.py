"""Certify an exact low-rank Gram completion near a rational seed."""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx, fmpq

root = Path(__file__).parent
prev = run_path(str(root/"ch39-local-memory-rank.py"))


def row_norm(matrix):
    return max(sum(abs(a) for a in row) for row in matrix.tolist())


def verify():
    ctx.prec = 192
    _, _, _, g0, g1, g2 = prev["MEM"]["coefficients"]()
    delta = g2-g1*g0.solve(g1)
    points = sorted(prev["LOCAL"]["shell_grid"](1)|{(17,17,17)})
    pairs = [(i,j) for i,j in combinations(range(19), 2) if prev["LOCAL"]["squared"](points[i],points[j])>56]
    larger_scale = {(11,11,25), (13,13,23), (17,17,17)}
    scales = [8192 if tuple(sorted(p)) in larger_scale else 2048 for p in points]
    seed = arb_mat([[arb(int(a))/10**17 for a in line.split()] for line in
                    (root/"ch39-local-gram-14.txt").read_text().splitlines()])
    n, rank = seed.nrows(), seed.ncols()
    assert (n, rank, len(pairs)) == (19,14,153)
    jac = arb_mat(len(pairs), n*rank)
    gram = seed*seed.transpose()
    error = arb(0)
    for t,(i,j) in enumerate(pairs):
        for k in range(rank):
            jac[t,i*rank+k] = seed[j,k]
            jac[t,j*rank+k] = seed[i,k]
        error = max(error, abs(gram[i,j]-arb(delta[i,j]*scales[i]*scales[j])))
    assert error < arb("1e-15")
    right = jac.transpose()*(jac*jac.transpose()).inv()
    bound = row_norm(right)
    assert bound < 3000
    radius = arb("1e-12")
    assert error + rank*(bound*radius)**2 < radius
    assert 2*rank*bound**2*radius < 1
    # Every exact solution produced in this box has the required contraction.
    enclosed = arb_mat([[(seed[i,k]+arb(0,bound*radius))/scales[i] for k in range(rank)] for i in range(n)])
    norm = row_norm(enclosed*enclosed.transpose())
    assert norm < arb("9e-7") and row_norm(g0.inv()) < 315000
    assert norm*315000 < (1-arb(fmpq(98,1000)))**2
    print("Gram rank at most",rank,"right inverse norm",bound.str(8),"residual",error.str(8))
    print("Row norm bound", norm.str(10),"auxiliary dimension at most",n+rank)


if __name__ == "__main__":
    verify()
