"""Four polynomial residues, the first graph band, and both relative edge readings."""
from pathlib import Path
from runpy import run_path

from sympy import Matrix, Poly, Rational, cancel, eye, ones, symbols, zeros

folder = Path(__file__).resolve().parent
cells = run_path(str(folder/"ch39-finite-shell-relative-index.py"))
graph = run_path(str(folder/"ch39-finite-shell-first-triplet.py"))
t, u, alpha = symbols("t u alpha")


def quotient(matrix, polynomial, variable):
    return matrix.applyfunc(lambda a: Poly(a.expand(), variable).rem(Poly(polynomial, variable)).as_expr())


def main():
    roots = (Rational(1), u, u**2)
    q3 = lambda m: quotient(m, u**2+u+1, u)
    # Each column gives the four simple-pole residues for one numerator 1,t,t².
    residues = Matrix([[0, 0, -1]]+[[r/3, r**2/3, Rational(1,3)] for r in roots])
    residues = q3(residues)
    assert q3(ones(1,4)*residues) == zeros(1,3)
    inv = Matrix([[0]+[r**k for r in roots] for k in (2,1,0)])
    assert q3(inv*residues) == eye(3)
    den = t**3-1
    # Pull back each of the three differentials under both rotations.
    for mobius, order in [(u*t, [0,2,3,1]), ((t+2)/(t-1), [1,0,3,2])]:
        cols = []
        for k in range(3):
            if mobius == u*t:
                numerator = u**(k+1)*t**k
            else:
                numerator = cancel(den*mobius**k*mobius.diff(t)/(mobius**3-1))
            pol = Poly(numerator, t)
            assert pol.degree() <= 2
            cols.append(Matrix([pol.nth(j) for j in range(3)]))
        assert q3(residues*Matrix.hstack(*cols)-residues.extract(order, range(3))) == zeros(4,3)
    points, edges, _ = cells["shell_cells"]()
    assert graph["points"] == points
    weight = dict(zip((0,1,3,4), (1,alpha,-alpha,-1)))
    V = Matrix([[weight[x] for x in p] for p in points])
    reduce5 = lambda m: quotient(m, alpha**2+alpha-1, alpha)
    L = graph["matrix"]({1:1, 2:1, 4:1})
    assert V*ones(4,1) == zeros(24,1)
    assert reduce5(L*V-(2-2*alpha)*V) == zeros(24,4)
    assert reduce5(V.T*V-16*(1+alpha**2)*(eye(4)-ones(4)/4)) == zeros(4)
    for label, sign in [(0,1),(4,-1)]:
        delta = []
        for a,b in edges:
            i,j = points[a].index(label), points[b].index(label)
            if i == j:
                continue
            row = eye(4)[j,:]-eye(4)[i,:]
            delta.append(row)
            assert V[b,:]-V[a,:] == sign*(1-alpha)*row
        B = Matrix.vstack(*delta)
        assert B.T*B == 8*eye(4)-2*ones(4)
    print("Residues intertwine both rotations; first-band gradients recover the two relative triplets with opposite sign.")


if __name__ == "__main__":
    main()
