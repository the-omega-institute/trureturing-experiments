"""Exact simplicial reflection and a fixed-volume witness on an untwisted double slab."""
from collections import Counter
from itertools import combinations, permutations, product

from flint import fmpq, fmpq_mat


def vertex(x):
    # Periods 4(e1+e2), 4e3, 4e4, 2(e1-e2).
    turns, h = divmod(x[0]-x[1],4)
    return (h,(x[1]+2*turns)%4,x[2]%4,x[3]%4)


def lift(v):
    h,b,c,d = v
    return (h+b,b,c,d)


def reflect(v):
    a,b,c,d = lift(v)
    return vertex((b,a,c,d))


def squared_gram(changed=None, shift=fmpq(1,64)):
    q = {}
    for i,j in combinations(range(5),2):
        d = j-i
        q[i,j] = q[j,i] = fmpq(d*(5-d),5) + (shift if (i,j) == changed else 0)
    for i in range(5):
        q[i,i] = 0
    g = fmpq_mat([[(q[0,i]+q[0,j]-q[i,j])/2 for i in range(1,5)]
                  for j in range(1,5)])
    for k in range(1,5):
        assert fmpq_mat([[g[i,j] for j in range(k)] for i in range(k)]).det() > 0
    assert g.det()/max(q.values())**4 > fmpq(1,16)
    return g.det()


def verify():
    cells = {}
    edges = set()
    incidence = Counter()
    changed = frozenset((vertex((0,0,0,0)),vertex((1,0,0,0))))
    for start in product(range(4), repeat=4):
        for order in permutations(range(4)):
            x = list(lift(start))
            path, heights = [start], [x[0]-x[1]]
            for i in order:
                x[i] += 1
                path.append(vertex(x))
                heights.append(x[0]-x[1])
            assert max(heights)-min(heights) == 1
            positive = min(heights)%4 in (0,1)
            key = frozenset(path)
            assert len(key) == 5 and key not in cells
            cells[key] = positive
            for i,j in combinations(range(5),2):
                edge = frozenset((path[i],path[j]))
                edges.add(edge)
                if edge == changed:
                    assert positive and j == i+1
                    incidence[i,j] += 1
    assert len(cells) == 6144 and sum(cells.values()) == 3072 and len(edges) == 3840
    for cell,positive in cells.items():
        assert cells[frozenset(map(reflect,cell))] != positive
        assert all(v == reflect(v) for v in cell if v[0] in (0,2))
    assert sum(incidence.values()) == 24
    det = squared_gram()
    assert det == fmpq(1,5) and det/fmpq(6,5)**4 == fmpq(125,1296)
    shift = fmpq(1,64)
    ratios = {squared_gram(pair,shift)/det for pair in incidence}
    assert ratios == {1+shift-3*shift**2/4}
    assert min(ratios) > 1
    print("Fixed seams, exchanged half-cell counts, edges:", 2, 3072, len(edges))
    print("Positive-half affected simplices and determinant ratio:", sum(incidence.values()), ratios)
    # In the two-choice model per half, F=-1 or +1 records the altered edge.
    def reflected_mean(weights):
        return sum((-1)**(i+j)*weights[i][j] for i,j in product(range(2),repeat=2)) / sum(map(sum,weights))
    one = fmpq(1)
    assert reflected_mean([[0,one],[one,0]]) == -1
    assert reflected_mean([[one/2,one],[one,one/2]]) == -fmpq(1,3)
    print("Rational four-geometry examples: fixed volume -1, soft volume -1/3.")


if __name__ == "__main__":
    verify()
