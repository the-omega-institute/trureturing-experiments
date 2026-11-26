"""精确验证线性层间跃迁的三循环选律、行列式和三个分离的局部根。"""
import sympy as sp


def verify():
    omega = (-1+sp.I*sp.sqrt(3))/2
    z, zbar = sp.symbols("z zbar")
    phase = sp.diag(1, omega, omega**2)
    variables = [sp.symbols(f"{letter}0:3") for letter in "abc"]
    constant, linear, conjugate = [sp.zeros(3) for _ in range(3)]
    terms = (constant, linear, conjugate)
    for i in range(3):
        for j in range(3):
            charge = (1+i-j) % 3
            terms[charge][i, j] = variables[charge][i]
    mixing = constant + z*linear + zbar*conjugate
    rotated = mixing.subs({z:omega*z, zbar:omega**2*zbar}, simultaneous=True)
    assert (rotated-omega*phase*mixing*phase.adjoint()).applyfunc(sp.simplify) == sp.zeros(3)
    det = sp.Poly(mixing.det(), z, zbar)
    assert set(det.monoms()) == {(0, 0), (3, 0), (0, 3), (1, 1)}
    assert det.coeff_monomial(z**3) == sp.prod(variables[1])
    assert det.coeff_monomial(1) == sp.prod(variables[0])

    radius = sp.Rational(1, 10**13)
    split = radius/10
    example = sp.Matrix([[z, 1, 0], [0, z, 1], [-split**3, 0, z]])
    f = sp.factor(example.det())
    for j in range(3):
        root = split*omega**j
        assert sp.simplify(f.subs(z, root)) == 0
        derivative = sp.diff(f, z).subs(z, root)
        assert sp.simplify(derivative*sp.conjugate(derivative)) == 9*split**4
    assert split < radius
    print("三循环实线性行列式的全部允许单项式", det.monoms())
    print("解析例的三个分离单根所在半径、包围半径", split, radius)


if __name__ == "__main__":
    verify()
