"""Exact weak-basis stabilizer tests for a three-generation rank-one contact."""

from sympy import I, Matrix, diag, eye, ones, symbols, zeros


def cyclic(grams, vector):
    """In dimension three a nonzero cyclic vector spans after two word extensions."""
    basis = vector
    for _ in range(2):
        extended = Matrix.hstack(basis, *(gram*basis for gram in grams))
        basis = Matrix.hstack(*extended.columnspace())
    return basis.rank() == 3


def hermitian_coordinates():
    real = symbols("a:9", real=True)
    a, b, c, d, e, f, g, h, k = real
    matrix = Matrix([[a, d+I*e, f+I*g], [d-I*e, b, h+I*k], [f-I*g, h-I*k, c]])
    return real, matrix


def linear_solutions(variables, template, constraints):
    entries = [part for constraint in constraints for entry in constraint
               for part in entry.expand(complex=True).as_real_imag()]
    equations = Matrix(entries).jacobian(variables)
    return [template.subs(dict(zip(variables, vector, strict=True))) for vector in equations.nullspace()]


def lepton_order(gram, weinberg, line):
    projector = line*line.H/(line.H*line)[0]
    assert gram == gram.H and weinberg == weinberg.T
    variables, hermitian = hermitian_coordinates()
    skew = I*hermitian
    tangent = linear_solutions(variables, skew, [
        skew*gram-gram*skew, skew*weinberg-weinberg*skew.conjugate(), skew*projector-projector*skew])
    if any((x.trace()-3*(projector*x).trace()).expand() != 0 for x in tangent):
        return 0
    transverse = linear_solutions(variables, hermitian, [
        hermitian*gram-gram*hermitian,
        hermitian*weinberg-weinberg*hermitian.conjugate(), hermitian*line])
    complement = eye(3)-projector
    if any((2*y-y.trace()*complement).applyfunc(lambda z: z.expand()) != zeros(3) for y in transverse):
        return 2
    return 1


def check():
    unit = eye(3)
    e1, e2, e3 = [unit[:, j] for j in range(3)]
    d = diag(1, 2, 3)
    dense = unit-2*ones(3)/3
    assert dense*dense.H == unit and all(element != 0 for element in dense)
    # Two exact triples of invertible Dirac Yukawas with the same dense supports.
    for yukawa, expected in ((dense, False), (dense*d, True)):
        assert yukawa.det() and all(element != 0 for element in yukawa)
        gram = yukawa*yukawa.H
        assert cyclic([gram, gram], e1) == expected
        assert lepton_order(gram, zeros(3), e1) == int(expected)
    assert cyclic([d, unit+ones(3)], e1)
    examples = [
        (d, zeros(3), e1+e2+e3, 1),
        (d, zeros(3), e1+e2, 0),
        (d, d, e1+e2+e3, 1),
        (d, d, e1, 2),
        (unit+e2*e3.H/2+e3*e2.H/2, d, e1, 1),
        (unit, unit, e1+I*e2, 0),
        (unit, unit, e1, 2),
        (d, diag(1, 2, 0), e1, 0),
        (d, diag(1, 2, 0), e1+e3, 2),
        (d, diag(1, 2, 0), e1+e2+e3, 1),
        (d, diag(1, 0, 0), e1+e2+e3, 1),
    ]
    change = diag(1, I, -1)*dense
    assert change.H*change == unit
    for gram, weinberg, line, expected in examples:
        assert lepton_order(gram, weinberg, line) == expected
        assert lepton_order(change*gram*change.H, change*weinberg*change.T,
                            (1+I)*change*line) == expected
    print("Dense Dirac textures and 11 exact lepton backgrounds verified in two weak bases.")


if __name__ == "__main__":
    check()
