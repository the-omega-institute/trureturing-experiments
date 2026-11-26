"""A finite signed-shell representative of unit weak-plane charge and its unreflected chart."""
from collections import Counter
from fractions import Fraction as F
from itertools import product


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def det(a, b, c):
    return sum(a[i]*(b[(i+1)%3]*c[(i+2)%3]-b[(i+2)%3]*c[(i+1)%3]) for i in range(3))


def verify():
    shell = [v for v in product(range(30), repeat=3) if dot(v, v) == 867]
    assert len(shell) == 19 and all(v[2] >= 0 for v in shell)
    axes = ((29, 1, 5), (5, 29, 1), (1, 5, 29))
    assert all(v in shell for v in axes)
    oriented = Counter()
    real_numerators = Counter()
    for signs in product((-1, 1), repeat=3):
        vertices = [tuple(s*x for x in v) for s, v in zip(signs, axes)]
        if det(*vertices) < 0:
            vertices[1], vertices[2] = vertices[2], vertices[1]
        a, b, c = vertices
        assert det(a, b, c) > 0
        real = 867+dot(a,b)+dot(b,c)+dot(c,a)
        assert real > 0
        real_numerators[real] += 1
        for i, v in enumerate(vertices):
            oriented[v, vertices[(i+1)%3]] += 1
    assert len(oriented) == 24
    assert all(count == 1 and oriented[b,a] == 1 for (a,b), count in oriented.items())
    assert real_numerators == {867+3*179: 2, 867-179: 6}
    # One winding in det(U(2)) has nontrivial endpoints in the simply connected factors.
    color = (F(-2, 3), F(1, 3), F(1, 3))
    weak = (F(1, 2), F(-1, 2))
    hyper = (F(-1, 3),)*3+(F(1, 2),)*2
    assert tuple(x+y for x,y in zip(color+weak, hyper)) == (-1,0,0,1,0)
    assert set(x%1 for x in color) == {F(1,3)} and set(x%1 for x in weak) == {F(1,2)}
    print('Signed-shell vertices, edges, oriented faces:', 6, len(oriented)//2, sum(real_numerators.values()))
    print('Trace-product real numerators and multiplicities:', real_numerators, 'positive determinant:', det(*axes))


if __name__ == '__main__':
    verify()
