"""Simple-cycle determinant derivative and spectral-top ordering on the Behrend graph."""
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
from sympy import Matrix, Poly, Rational, cancel, chebyshevt, symbols

shell = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))


def check_polynomials():
    t = symbols("t")
    for k in range(1, 7):
        polynomial = chebyshevt(k, t)
        for value in (polynomial, polynomial.diff(t)):
            coefficients = Poly(value.subs(t, t+Rational(99,100)), t).all_coeffs()
            assert min(coefficients) > 0


def check_expansions():
    # The two triangles of the four-site diamond have a charged outer square.
    # These unequal rational gains and diagonal potentials also test cofactor signs.
    w = symbols("w", nonzero=True)
    z = symbols("z")
    m = Matrix.diag(2, -3, 5, 1)
    for i,j,a,k in [(0,1,2,1), (1,2,7,1), (2,3,3,0), (3,0,11,0), (1,3,5,0)]:
        m[i,j], m[j,i] = a*w**k, a*w**(-k)
    pencil = z*Matrix.eye(4)-m
    vertices = tuple(range(4))
    def complement(path):
        rem = [i for i in vertices if i not in path]
        return pencil.extract(rem, rem).det()
    def weight(path):
        result = 1
        for i,j in zip(path, path[1:]):
            result *= m[i,j]
        return result
    derivative = 0
    for size in (3,4):
        for group in combinations(vertices, size):
            for order in permutations(group[1:]):
                cycle = (group[0], *order, group[0])
                term = weight(cycle)
                if term != 0:
                    derivative -= term.diff(w)*complement(group)
    assert cancel(derivative-pencil.det().diff(w)) == 0
    adj = pencil.adjugate()
    for i,j in permutations(vertices, 2):
        paths = 0
        for k in range(3):
            for middle in permutations((v for v in vertices if v not in (i,j)), k):
                path = (i, *middle, j)
                paths += weight(path)*complement(path)
        assert cancel(paths-adj[i,j]) == 0


def verify():
    ctx.prec = 192
    check_polynomials()
    check_expansions()
    nodes, edges = shell["geometry"]()
    adjacency = [set() for p in nodes]
    charge = {}
    for i,j,_,k in edges:
        adjacency[i].add(j)
        adjacency[j].add(i)
        charge[i,j], charge[j,i] = k, -k
    visited, pending = {0}, [0]
    while pending:
        i = pending.pop()
        for j in adjacency[i]-visited:
            visited.add(j)
            pending.append(j)
    assert len(visited) == len(nodes) and sum(abs(k) for *_,k in edges) == 6
    triangle = [nodes.index(p) for p in [(21,51,69), (33,33,75), (39,39,69), (21,51,69)]]
    assert sum(charge[i,j] for i,j in zip(triangle,triangle[1:])) == 1
    assert acb.zeta_zero(1).imag > 14
    assert (4*14**2-1)*100 > 99*(4*14**2+1)
    maxima, determinants = [], []
    for n in [1, *range(11485,11492)]:
        rho = acb.zeta_zero(n)
        matrix = shell["matrices"](nodes, edges, (rho-1)/rho)
        eig = matrix.eig(multiple=True)
        top = max(v.real for v in eig)
        assert len(eig) == len(nodes) and sum(not v.real < top for v in eig) == 1
        pencil = -matrix
        for i in range(len(nodes)):
            pencil[i,i] += 1
        det = pencil.det()
        assert det.imag.contains(0) and det.real > 0
        maxima.append(top)
        determinants.append(det.real)
    assert all(a < b for a,b in zip(maxima, maxima[1:]))
    assert all(a > b for a,b in zip(determinants, determinants[1:]))
    gaps = [b-a for a,b in zip(maxima[1:], maxima[2:])]
    print("First and last high-zero spectral tops", maxima[1].str(29), maxima[-1].str(29))
    print("Smallest adjacent increase", min(gaps).str(13))
    print("Determinant at first and last high zero", determinants[1].str(27), determinants[-1].str(27))


if __name__ == "__main__":
    verify()
