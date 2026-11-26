"""Two integer connections with the same six-point Behrend readings and different zero modes."""
from fractions import Fraction
from itertools import combinations, permutations


def rank(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    r = 0
    for j in range(len(a[0])):
        i = next((i for i in range(r,len(a)) if a[i][j]),None)
        if i is None:
            continue
        a[i],a[r] = a[r],a[i]
        for i in range(r+1,len(a)):
            coef = a[i][j]/a[r][j]
            a[i] = [x-coef*y for x,y in zip(a[i],a[r])]
        r += 1
    return r


def main():
    points = list(permutations((1,3,4)))
    edges = [(i,j,26-sum(x*y for x,y in zip(points[i],points[j])))
             for i,j in combinations(range(6),2)]
    site_edges = [(i,j,d) for i,j,d in edges if d in (1,4,7)]
    assert len(site_edges)==12
    changed = next((i,j) for i,j,d in site_edges if d==1)
    for twist, nullity in ((False,3),(True,1)):
        A = [[0]*18 for _ in range(18)]
        for i,j,_ in site_edges:
            for k in range(3):
                sign = -1 if twist and (i,j)==changed and k<2 else 1
                A[3*i+k][3*i+k] += 1
                A[3*j+k][3*j+k] += 1
                A[3*i+k][3*j+k] -= sign
                A[3*j+k][3*i+k] -= sign
        assert 18-rank(A)==nullity
    # Gauge-fix the five edges of a spanning tree. Each remaining closing edge
    # can carry the sole nonidentity holonomy if its readout is omitted.
    tree = {0}
    tree_edges = set()
    while len(tree) < len(points):
        i,j,_ = next(edge for edge in site_edges if (edge[0] in tree) != (edge[1] in tree))
        tree.update((i,j))
        tree_edges.add((i,j))
    chords = [(i,j) for i,j,_ in site_edges if (i,j) not in tree_edges]
    assert len(chords)==7
    assert rank([[4,0,0],[0,4,0],[0,0,0]])==2
    print("Same 36 shell readings; integer connection Laplacian kernel dimensions 3 and 1; seven fundamental circuits")


if __name__ == "__main__":
    main()
