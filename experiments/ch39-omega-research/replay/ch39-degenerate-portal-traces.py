"""Exact finite-word closure and singular spectral padding for degenerate portals."""

from sympy import I, Matrix, Rational, diag, eye, interpolate, symbols, zeros


def word_closure(generators):
    columns = []
    pivots = []

    def retain(word):
        vector = list(word)
        for pivot, basis in pivots:
            scale = vector[pivot]
            if scale:
                vector = [a-scale*b for a, b in zip(vector, basis)]
        pivot = next((i for i, value in enumerate(vector) if value), None)
        if pivot is None:
            return
        scale = vector[pivot]
        pivots.append((pivot, [value/scale for value in vector]))
        columns.append(word)

    retain(eye(generators[0].rows))
    index = 0
    while index < len(columns):
        word = columns[index]
        for generator in generators:
            retain(word*generator)
        index += 1
    return columns


def check():
    portal = Matrix([[0, 1, -I], [1, 0, 1], [I, 1, 0]])
    projectors = [diag(*[int(i == j) for i in range(3)]) for j in range(3)]
    assert len(word_closure([*projectors, portal])) == 9
    word = projectors[0]*portal*projectors[1]*portal*projectors[2]*portal
    assert word.trace() == I
    assert word.trace() + word.trace().conjugate() == 0
    assert portal.eigenvals().get(0) == 1
    assert portal**3 == 3*portal

    # Distinct nonzero spectra give an invertible padded moment system.
    values = [-3, -1, 2, 4]
    system = Matrix([[v**k for v in values] for k in range(3, 3+len(values))])
    assert system.det() != 0
    x = symbols("x")
    q = interpolate([(v, Rational(1, v*v)) for v in values], x)
    assert all(v**3*q.subs(x, v) == v for v in [0, *values])

    # The six physical coordinates occur in three mass spaces of rank two.
    P = [diag(p, p) for p in projectors]
    left, right = diag(portal, portal), diag(portal, portal.conjugate())
    assert sum(P, zeros(6)) == eye(6)
    closed = [((P[0]*f)*(P[1]*f)*(P[2]*f)).trace() for f in (left, right)]
    assert closed == [2*I, 0]
    # Conjugation doubling gives two copies of each oriented triangle.
    # Reorder the three-dimensional summands to match the two doubled tuples.
    ordering = [*range(3), *range(6, 12), *range(3, 6)]
    exchange = eye(12)[ordering, :]
    for generator_left, generator_right in [*zip(P, P), (left, right)]:
        doubled_left = diag(generator_left, generator_left.conjugate())
        doubled_right = diag(generator_right, generator_right.conjugate())
        assert exchange*doubled_left*exchange.T == doubled_right
    print("Simple algebra dimension 9; rank-two mass-space triangle traces", closed)


if __name__ == "__main__":
    check()
