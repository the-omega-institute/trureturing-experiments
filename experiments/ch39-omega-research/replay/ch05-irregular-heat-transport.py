"""Exact generalized-moment collisions for fixed irregular cutoff schedules."""

from math import lcm, prod
from pathlib import Path
from runpy import run_path

from sympy import Matrix, Rational


helpers = run_path(str(Path(__file__).with_name("ch05-heat-trace-transport.py")))
moment, transport = helpers["moment"], helpers["transport"]


def tail_example(degrees, endpoint):
    n = len(degrees) - int(endpoint)
    assert n > 0 and n % 2 == 0
    denominator = lcm(*(a.q for a in degrees))
    nodes = [Rational(j+1, n+2)**denominator for j in range(n+1)]
    powers = [0, *degrees[:n-1]]
    values = Matrix([[x**a for x in nodes] for a in powers])
    minors = [values[:, [j for j in range(n+1) if j != i]].det() for i in range(n+1)]
    assert all(w > 0 for w in minors) and sum(minors[::2]) == sum(minors[1::2])
    z = sum(minors[::2])
    base = [{x: w/z for x, w in zip(nodes[p::2], minors[p::2])} for p in (1, 0)]
    assert all(moment(base[0], a) == moment(base[1], a) for a in powers)
    assert moment(base[0], degrees[n-1]) != moment(base[1], degrees[n-1])
    distance = transport(*base)
    assert distance > 0
    for level in (2, 4, 8):
        r = Rational(1, level)**denominator
        pair = [{r*x: p/(1+int(endpoint)) for x, p in row.items()} for row in base]
        if endpoint:
            for row in pair:
                row[1] = Rational(1, 2)
        assert len(pair[0]) == (len(degrees)+int(endpoint))//2
        assert transport(*pair) == r*distance/(1+int(endpoint))
        for a in degrees:
            error = moment(pair[0], a) - moment(pair[1], a)
            assert error == r**a*(moment(base[0], a)-moment(base[1], a))/(1+int(endpoint))


def check():
    for exponents, endpoint in [((1, 2, 5, 8), False), ((1, 2, 10), True),
                                ((1, 4, 10), True), ((Rational(1, 2), Rational(7, 2)), False),
                                ((Rational(1, 2), Rational(7, 2), 9), True)]:
        tail_example(tuple(map(Rational, exponents)), endpoint)
    # For the logarithmic derivative basis, check the row-difference recursion
    # at integer exponents using exact integrals on the original x scale.
    nodes = [Rational(1, 8), Rational(1, 4), Rational(3, 4)]
    beta = [Rational(-1, 2), Rational(3, 2), Rational(7, 2)]
    matrix = Matrix([[x**b for b in beta] for x in nodes])
    differences = [b-beta[0] for b in beta[1:]]
    integrated = Matrix([[(v**b-u**b)/b for b in differences] for u, v in zip(nodes, nodes[1:])])
    assert (matrix.det()-prod(x**beta[0] for x in nodes)*prod(differences)*integrated.det()).simplify() == 0
    print("Five irregular schedules: exact signed minors, error scaling, transport and determinant recursion verified.")


if __name__ == "__main__":
    check()
